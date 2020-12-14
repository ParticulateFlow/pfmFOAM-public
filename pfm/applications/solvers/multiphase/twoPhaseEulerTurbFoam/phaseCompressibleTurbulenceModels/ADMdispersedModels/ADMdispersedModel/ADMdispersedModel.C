/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "ADMdispersedModel.H"
#include "mathematicalConstants.H"
#include "twoPhaseSystem.H"
#include "fvOptions.H"
#include "fixedValueFvsPatchFields.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RASModels::ADMdispersedModel::ADMdispersedModel
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& phase,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity
    <
        RASModel<EddyDiffusivity<phaseCompressibleTurbulenceModel>>
    >
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        phase,
        propertiesName
    ),

    phase_(phase),

    regularizationModel_
    (
        ADMdispersedModels::regularizationModel::New
        (
            coeffDict_,
            alpha_
        )
    ),
    frictionalStressModel_
    (
        ADMdispersedModels::frictionalStressModel::New
        (
            coeffDict_
        )
    ),

    alphaMax_("alphaMax", dimless, coeffDict_),
    alphaMinFriction_
    (
        "alphaMinFriction",
        dimless,
        coeffDict_
    ),
    residualAlpha_
    (
        "residualAlpha",
        dimless,
        coeffDict_
    ),

    maxK_
    (
        "maxK",
        dimensionSet(0,2,-2,0,0),
        coeffDict_.lookupOrDefault<scalar>("maxK",25.0)
    ),

    deconOrder_
    (
        "deconOrder",
        dimensionSet(0,0,0,0,0),
        coeffDict_.lookupOrDefault<scalar>("deconOrder", 7)
    ),

    maxNut_
    (
        "maxNut",
        dimensionSet(0,2,-1,0,0),
        coeffDict_.lookupOrDefault<scalar>("maxNut",1000)
    ),

    R1ADM_
    (
        IOobject
        (
            IOobject::groupName("RADM", phase.name()),
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),

    Dstar_
    (
        IOobject
        (
            IOobject::groupName("Dstar", phase.name()),
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
        dimensionedSymmTensor("zero", dimensionSet(0, 0, -1, 0, 0),
                          symmTensor(0.0,0.0,0.0,0.0,0.0,0.0)),
        zeroGradientFvPatchField<scalar>::typeName
    ),

    alpha1star_
    (
        IOobject
        (
            "alpha1star",
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("value", dimensionSet(0, 0, 0, 0, 0), 0.5),
        zeroGradientFvPatchField<scalar>::typeName
    ),

    alphaP2Mean_
    (
        IOobject
        (
            "alphaP2Mean",
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.0),
        zeroGradientFvPatchField<scalar>::typeName
    ),


    U1star_
    (
        IOobject
        (
            "U1star",
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),

    nuFric_
    (
        IOobject
        (
            IOobject::groupName("nuFric", phase.name()),
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("zero", dimensionSet(0, 2, -1, 0, 0), 0.0),
        zeroGradientFvPatchField<scalar>::typeName
    ),

    pf_
    (
        IOobject
        (
            IOobject::groupName("pf", phase.name()),
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("zero", dimensionSet(1, -1, -2, 0, 0), 0.0),
        zeroGradientFvPatchField<scalar>::typeName
     ),

    k_
    (
        IOobject
        (
            IOobject::groupName("k", phase.name()),
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("zero", dimensionSet(0, 2, -2, 0, 0), 0.0)
     ),

    filterPtr_(LESfilter::New(U.mesh(), coeffDict_)),
    filter_(filterPtr_())

{
    if (type == typeName)
    {
        printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RASModels::ADMdispersedModel::~ADMdispersedModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::RASModels::ADMdispersedModel::read()
{
    if
    (
        eddyViscosity
        <
            RASModel<EddyDiffusivity<phaseCompressibleTurbulenceModel>>
        >::read()
    )
    {
        alphaMax_.readIfPresent(coeffDict());
        alphaMinFriction_.readIfPresent(coeffDict());
        residualAlpha_.readIfPresent(coeffDict());
        maxK_.readIfPresent(coeffDict());
        deconOrder_.readIfPresent(coeffDict());
        regularizationModel_->read();
        frictionalStressModel_->read();

        return true;
    }
    else
    {
        return false;
    }
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::ADMdispersedModel::k() const
{
    return k_;
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::ADMdispersedModel::epsilon() const
{
    NotImplemented;
    return nut_;
}


Foam::tmp<Foam::volSymmTensorField>
Foam::RASModels::ADMdispersedModel::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                IOobject::groupName("R", U_.group()),
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
          - (nut_)*dev(twoSymm(fvc::grad(U_)))
          + symm(R1ADM_)
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::ADMdispersedModel::pPrime() const
{
    const volScalarField& rho = phase_.rho();
    tmp<volScalarField> tda(phase_.d());
    const volScalarField& da = tda();
    
    tmp<volScalarField> tpPrime
    (
        filter_(
            frictionalStressModel_->frictionalPressurePrime
            (
                phase_,
                alpha1star_,
                alphaMinFriction_,
                alphaMax_,
                da,
                rho,
                dev(Dstar_)
            )
        )*pos(alpha_-alphaMinFriction_)
    );
    
    volScalarField::Boundary& bpPrime =
        tpPrime.ref().boundaryFieldRef();

    forAll(bpPrime, patchi)
    {
        if (!bpPrime[patchi].coupled())
        {
            bpPrime[patchi] = 0;
        }
    }

    return tpPrime;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::RASModels::ADMdispersedModel::pPrimef() const
{
    return fvc::interpolate(pPrime());
}

Foam::tmp<Foam::volScalarField>
Foam::RASModels::ADMdispersedModel::normalStress() const
{
    const volScalarField& rho = phase_.rho();
    tmp<volScalarField> tda(phase_.d());
    const volScalarField& da = tda();

    tmp<volScalarField> tNormalStress
    (
       pos(alpha_ - alphaMinFriction_)
      *filter_
       (
           frictionalStressModel_->frictionalPressure
           (
               phase_,
               alpha1star_,
               alphaMinFriction_,
               alphaMax_,
               da,
               rho,
               dev(Dstar_)
           )
       )
    );

    return tNormalStress;
}


Foam::tmp<Foam::volSymmTensorField>
Foam::RASModels::ADMdispersedModel::devRhoReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                IOobject::groupName("devRhoReff", U_.group()),
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           // Frictional stress
          - (rho_*nut_)
           *dev(twoSymm(fvc::grad(U_)))
           // Reynolds stress
          + symm(R1ADM_)
        )
    );
}


Foam::tmp<Foam::fvVectorMatrix>
Foam::RASModels::ADMdispersedModel::divDevRhoReff
(
    volVectorField& U
) const
{
    return
    (
      - fvm::laplacian(rho_*nut_, U)
      - fvc::div
        (
            // frictional stress
            (rho_*nut_)*dev2(T(fvc::grad(U)))
        )
     - pos(alpha_ - residualAlpha_)*fvc::div
        (// Reynolds stress
          - symm(R1ADM_)
        )
      //regularization
      + regularizationModel_->regTerm(alpha_,rho_,k_,U,U1star_)
    );
}

void Foam::RASModels::ADMdispersedModel::boundNormalStress
(
    volTensorField& R
) const
{
    scalar kMin = 1.e-7;
    scalar kMax = maxK_.value();

    R.max
    (
        dimensionedTensor
        (
            "zero",
            R.dimensions(),
            tensor
            (
                  kMin, - kMax, - kMax,
                - kMax,   kMin, - kMax,
                - kMax, - kMax,   kMin
            )
        )
    );
    
    R.min
    (
        dimensionedTensor
        (
            "zero",
            R.dimensions(),
            tensor
            (
                  kMax, kMax, kMax,
                  kMax, kMax, kMax,
                  kMax, kMax, kMax
            )
        )
    );
}


void Foam::RASModels::ADMdispersedModel::correct()
{
    // Local references
    volScalarField alpha(max(alpha_, scalar(0)));
    const volScalarField& rho = phase_.rho();
    const volVectorField& U = U_;

    tmp<volScalarField> tda(phase_.d());
    const volScalarField& da = tda();
    
    // ADM
    alpha1star_  = alpha;
    volScalarField alpha1starT = alpha;
    U1star_      = U;
    
    
    for (int i = 0; i < int(deconOrder_.value()); i++) {
        alpha1starT   = filter_(alpha1star_);
        alpha1starT.min(alphaMax_.value());
        alpha1starT.max(1.0e-6);
        U1star_      += U - filter_(alpha1star_*U1star_)/alpha1starT;
        alpha1star_  += alpha - alpha1starT;
        alpha1star_.min(alphaMax_.value());
        alpha1star_.max(1.0e-6);
    }
    alpha1star_.correctBoundaryConditions();
    U1star_.correctBoundaryConditions();
    
    volScalarField a1sF = filter_(alpha1star_);
    a1sF.min(alphaMax_.value());
    a1sF.max(1.0e-6);
    volVectorField U1sF = filter_(alpha1star_*U1star_)/a1sF;
    
    // Compute alphaP2Mean for virutal mass model of de Wilde
    alphaP2Mean_ = filter_(sqr(alpha1star_)) - sqr(a1sF);
    alphaP2Mean_.max(sqr(1.0e-6));
    
    // Compute Reynolds stress tensor
    R1ADM_ = rho*(
                    filter_(alpha1star_ * U1star_ * U1star_)
                  - a1sF * U1sF * U1sF
                );
    
    // limit Reynolds stress
    boundNormalStress(R1ADM_);
    R1ADM_.correctBoundaryConditions();
    
    // compute turbulent kinetic energy
    k_ = tr(R1ADM_)/(rho*max(a1sF,residualAlpha_));

    // compute derivative of Ustar for frictional model
    tmp<volTensorField> tgradU(fvc::grad(U1star_));
    const volTensorField& gradU(tgradU());
    Dstar_ = symm(gradU);

    // Frictional pressure
    pf_ = frictionalStressModel_->frictionalPressure
    (
        phase_,
        alpha1star_,
        alphaMinFriction_,
        alphaMax_,
        da,
        rho,
        dev(Dstar_)
    );
    
    nuFric_ = frictionalStressModel_->nu
    (
        phase_,
        alpha1star_,
        alphaMinFriction_,
        alphaMax_,
        pf_/rho,
        da,
        dev(Dstar_)
    );
    // filter frictional pressure and viscosity
    pf_ = filter_(pf_);
    nut_ = filter_(nuFric_)*pos(alpha - alphaMinFriction_);
    // Limit viscosity and add frictional viscosity
    nut_.min(maxNut_.value());
    nut_.max(0.);
    
    Info<< "ADM (dispersed):" << nl
        << "    max(nut) = " << max(nut_).value() << nl
        << "    max(k)   = " << max(k_).value()   << endl;

    if (debug)
    {
        Info<< typeName << ':' << nl
            << "    max(nut) = " << max(nut_).value() << endl;
    }
}


// ************************************************************************* //

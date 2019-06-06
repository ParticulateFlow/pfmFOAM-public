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

#include "ADMcontinuousModel.H"
#include "mathematicalConstants.H"
#include "twoPhaseSystem.H"
#include "fvOptions.H"
#include "fixedValueFvsPatchFields.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RASModels::ADMcontinuousModel::ADMcontinuousModel
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
        ADMcontinuousModels::regularizationModel::New
        (
            coeffDict_,
            alpha_
        )
    ),

    alphaMax_("alphaMax", dimless, coeffDict_),

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

    Cmu_
    (
        "Cmu",
        dimensionSet(0,0,0,0,0),
        coeffDict_.lookupOrDefault<scalar>("Cmu",0.4)
    ),

    R2ADM_
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

    alpha2star_
    (
        IOobject
        (
            "alpha2star",
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("value", dimensionSet(0, 0, 0, 0, 0), 0.5),
        // Set Boundary condition
        zeroGradientFvPatchField<scalar>::typeName
    ),

    U2star_
    (
        IOobject
        (
            "U2star",
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
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

Foam::RASModels::ADMcontinuousModel::~ADMcontinuousModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::RASModels::ADMcontinuousModel::read()
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
        residualAlpha_.readIfPresent(coeffDict());
        maxK_.readIfPresent(coeffDict());
        deconOrder_.readIfPresent(coeffDict());
        Cmu_.readIfPresent(coeffDict());
        regularizationModel_->read();

        return true;
    }
    else
    {
        return false;
    }
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::ADMcontinuousModel::k() const
{
    return k_;
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::ADMcontinuousModel::epsilon() const
{
    NotImplemented;
    return nut_;
}


Foam::tmp<Foam::volSymmTensorField>
Foam::RASModels::ADMcontinuousModel::R() const
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
            symm(R2ADM_)
        )
    );
}

Foam::tmp<Foam::volSymmTensorField>
Foam::RASModels::ADMcontinuousModel::devRhoReff() const
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
           // Reynolds stress
            symm(R2ADM_)
        )
    );
}


Foam::tmp<Foam::fvVectorMatrix>
Foam::RASModels::ADMcontinuousModel::divDevRhoReff
(
    volVectorField& U
) const
{
    return
    (
      - fvm::laplacian(rho_*nut_*scalar(0), U)
      - fvc::div
        (
            // frictional stress
            (rho_*nut_*scalar(0))*dev2(T(fvc::grad(U)))
            // Reynolds stress
          - symm(R2ADM_)
        )
      //regularization
      + regularizationModel_->regTerm(alpha_,rho_,k_,U,U2star_)
    );
}

void Foam::RASModels::ADMcontinuousModel::boundNormalStress
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

void Foam::RASModels::ADMcontinuousModel::correct()
{
    // Local references
    volScalarField alpha(max(alpha_, scalar(0)));
    const volScalarField& rho = phase_.rho();
    const volVectorField& U = U_;
    
    volScalarField cellVolume
    (
        IOobject
        (
            "cellVolume",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("one", dimLength*dimLength*dimLength, 1)
    );
    // gradient of velocity
    tmp<volTensorField> tgradU(fvc::grad(U_));
    const volTensorField& gradU(tgradU());
    volSymmTensorField D(dev(symm(gradU)));
    
    // ADM
    volScalarField alpha1      = scalar(1.0) - alpha;
    volScalarField alpha1star  = alpha1;
    volScalarField alpha1starT = alpha1;
    volScalarField alpha2starT = alpha;
    alpha2star_  = alpha;
    U2star_      = U;
    
    
    for (int i = 0; i < int(deconOrder_.value()); i++) {
        alpha1starT   = filter_(alpha1star);
        alpha1starT.min(alphaMax_.value());
        alpha1starT.max(1.0e-6);
        alpha2starT   = scalar(1.0)  - alpha1starT;
        U2star_      += U - filter_(alpha2star_*U2star_)/alpha2starT;
        alpha1star   += alpha1 - alpha1starT;
        alpha1star.min(alphaMax_.value());
        alpha1star.max(1.0e-6);
        alpha2star_  = scalar(1.0) - alpha1star;
    }
    alpha2star_.correctBoundaryConditions();
    U2star_.correctBoundaryConditions();
    
    volScalarField a2sF = filter_(alpha2star_);
    a2sF.min(1.0);
    a2sF.max(1.0 - alphaMax_.value());
    volVectorField U2sF = filter_(alpha2star_*U2star_)/a2sF;
    
    // Compute Reynolds stress tensor
    R2ADM_ = rho
               *(
                    filter_(alpha2star_ * U2star_ * U2star_)
                  - a2sF * U2sF * U2sF
                );
    // limit Reynolds stress
    boundNormalStress(R2ADM_);
    R2ADM_.correctBoundaryConditions();
    
    // compute turbulent kinetic energy
    k_ = tr(R2ADM_)/(rho*alpha);
    
    cellVolume.ref() = mesh_.V();
    // regularized nut
    nut_ =  Cmu_*alpha*pow(cellVolume,1.0/3.0)
          * (sqrt(k_) + pow(cellVolume,1.0/3.0)*sqrt(D&&D));
    nut_.correctBoundaryConditions();
    
    Info<< "ADM (continuous):" << nl
        << "    max(nut) = " << max(nut_).value() << nl
        << "    max(k)   = " << max(k_).value()   << endl;
    
    if (debug)
    {
        Info<< typeName << ':' << nl
            << "    max(nut) = " << max(nut_).value() << endl;
    }
}


// ************************************************************************* //

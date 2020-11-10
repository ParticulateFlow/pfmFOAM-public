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

#include "SATFMdispersedModel.H"
#include "mathematicalConstants.H"
#include "twoPhaseSystem.H"
#include "wallDist.H"
#include "uniformDimensionedFields.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RASModels::SATFMdispersedModel::SATFMdispersedModel
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
    frictionalStressModel_
    (
        kineticTheoryModels::frictionalStressModel::New
        (
            coeffDict_
        )
    ),

    equilibriumK_(coeffDict_.lookup("equilibriumK")),
    equilibriumPhiP2_(coeffDict_.lookup("equilibriumPhiP2")),
    dynamicAdjustment_(coeffDict_.lookup("dynamicAdjustment")),
    anIsoTropicNut_(coeffDict_.lookup("anIsoTropicNut")),
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

    maxNut_
    (
        "maxNut",
        dimensionSet(0,2,-1,0,0),
        coeffDict_.lookupOrDefault<scalar>("maxNut",10)
    ),

    xiPhiSolidScalar_
    (
        "xiPhiSolidScalar",
        dimensionSet(0,0,0,0,0),
        coeffDict_.lookupOrDefault<scalar>("xiPhiS",-0.1)
    ),

    xiPhiDivUScalar_
    (
        "xiPhiDivUScalar",
        dimensionSet(0,0,0,0,0),
        coeffDict_.lookupOrDefault<scalar>("xiPhiDivU",-0.15)
    ),

    xiGSScalar_
    (
        "xiGSScalar",
        dimensionSet(0,0,0,0,0),
        coeffDict_.lookupOrDefault<scalar>("xiGS",0.9)
    ),

    CmuScalar_
    (
        "CmuScalar",
        dimensionSet(0,0,0,0,0),
        coeffDict_.lookupOrDefault<scalar>("Cmu",0.2)
    ),

    CmuWScalar_
    (
        "CmuWScalar",
        dimensionSet(0,0,0,0,0),
        coeffDict_.lookupOrDefault<scalar>("CmuW",0.3)
    ),

    CphiSscalar_
    (
        "CphiSscalar",
        dimensionSet(0,0,0,0,0),
        coeffDict_.lookupOrDefault<scalar>("CphiS",0.3)
    ),
    CepsScalar_
    (
        "CepsScalar",
        dimensionSet(0,0,0,0,0),
        coeffDict_.lookupOrDefault<scalar>("Ceps",1.0)
    ),

    CpScalar_
    (
        "CpScalar",
        dimensionSet(0,0,0,0,0),
        coeffDict_.lookupOrDefault<scalar>("Cp",0.4)
    ),

    sigma_
    (
        "sigmaD",
        dimensionSet(0,0,0,0,0),
        coeffDict_.lookupOrDefault<scalar>("Sigma",2.0)
    ),

    gN_
    (
        "g",
        dimensionSet(0,1,-2,0,0),
        coeffDict_.lookupOrDefault<vector>("g",vector(0,0,-9.81))
    ),

    maxK_
    (
        "maxK",
        dimensionSet(0,2,-2,0,0),
        coeffDict_.lookupOrDefault<scalar>("maxK",25.0)
    ),

    ut_
    (
        "ut",
        dimensionSet(0,1,-1,0,0),
        coeffDict_.lookupOrDefault<scalar>("ut",1.0)
    ),

    k_
    (
        IOobject
        (
            IOobject::groupName("k", phase.name()),
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
        dimensionedScalar("zero", dimensionSet(0, 2, -1, 0, 0), 0.0)
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
        dimensionedScalar("zero", dimensionSet(1, -1, -2, 0, 0), 0.0)
    ),

    xiPhiS_
    (
        IOobject
        (
            "xiPhiS",
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedVector("value", dimensionSet(0, 0, 0, 0, 0), vector(-0.1,-0.1,-0.1)),
        zeroGradientFvPatchField<vector>::typeName
    ),

    xiPhiDivU_
    (
        IOobject
        (
            "xiPhiDivU",
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("value", dimensionSet(0, 0, 0, 0, 0), -0.15),
        zeroGradientFvPatchField<scalar>::typeName
    ),

   xiPhi2DivU_
    (
        IOobject
        (
            "xiPhi2DivU",
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("value", dimensionSet(0, 0, 0, 0, 0), -0.15),
        zeroGradientFvPatchField<scalar>::typeName
    ),

    xiUU_
    (
        IOobject
        (
            "xiUU",
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
        dimensionedTensor("value", dimensionSet(0, 0, 0, 0, 0), tensor(1,0,0, 0,1,0, 0,0,1)),
        zeroGradientFvPatchField<tensor>::typeName
    ),


    xiPhiGG_
    (
        IOobject
        (
            "xiPhiGG",
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("value", dimensionSet(0, 0, 0, 0, 0), 1.0),
        zeroGradientFvPatchField<scalar>::typeName
    ),

    xiGS_
    (
        IOobject
        (
            "xiGS",
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("value", dimensionSet(0, 0, 0, 0, 0), 0.9),
        zeroGradientFvPatchField<scalar>::typeName
    ),

    xiGatS_
    (
        IOobject
        (
            "xiGatS",
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("value", dimensionSet(0, 0, 0, 0, 0), 1.0),
        zeroGradientFvPatchField<scalar>::typeName
    ),

    alphaP2Mean_
    (
        IOobject
        (
            IOobject::groupName("alphaP2Mean", phase.name()),
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),

    Cmu_
    (
        IOobject
        (
            IOobject::groupName("Cmu", phase.name()),
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("value", dimensionSet(0, 0, 0, 0, 0), 1.0e-2)
    ),

    Ceps_
    (
        IOobject
        (
            IOobject::groupName("Ceps", phase.name()),
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("value", dimensionSet(0, 0, 0, 0, 0), 1.0)
    ),

    Cp_
    (
        IOobject
        (
            IOobject::groupName("Cp", phase.name()),
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("value", dimensionSet(0, 0, 0, 0, 0), 1.0)
    ),

    CphiS_
    (
        IOobject
        (
            IOobject::groupName("CphiS", phase.name()),
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("value", dimensionSet(0, 0, 0, 0, 0), 0.1)
    ),

    deltaF_
    (
        IOobject
        (
            IOobject::groupName("deltaF", phase.name()),
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("value", dimensionSet(0, 1, 0, 0, 0), 1.e-2)
    ),

    lm_
    (
        IOobject
        (
            IOobject::groupName("lm", phase.name()),
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("value", dimensionSet(0, 1, 0, 0, 0), 1.e-2)
    ),

    R1_
    (
        IOobject
        (
            IOobject::groupName("R", phase.name()),
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),

    shearProd_
    (
        IOobject
        (
            IOobject::groupName("shearProd", phase.name()),
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
        dimensionedVector("zero", dimensionSet(0, 2, -3, 0, 0),
                           vector(0.0,0.0,0.0)),
        zeroGradientFvPatchField<vector>::typeName
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

Foam::RASModels::SATFMdispersedModel::~SATFMdispersedModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::RASModels::SATFMdispersedModel::read()
{
    if
    (
        eddyViscosity
        <
            RASModel<EddyDiffusivity<phaseCompressibleTurbulenceModel>>
        >::read()
    )
    {
        coeffDict().lookup("equilibriumK") >> equilibriumK_;
        coeffDict().lookup("equilibriumPhiP2") >> equilibriumPhiP2_;
        coeffDict().lookup("dynamicAdjustment") >> dynamicAdjustment_;
        coeffDict().lookup("anIsoTropicNut") >> anIsoTropicNut_;
        alphaMax_.readIfPresent(coeffDict());
        alphaMinFriction_.readIfPresent(coeffDict());
        xiPhiSolidScalar_.readIfPresent(coeffDict());
        xiPhiDivUScalar_.readIfPresent(coeffDict());
        xiGSScalar_.readIfPresent(coeffDict());
        CmuScalar_.readIfPresent(coeffDict());
        CmuWScalar_.readIfPresent(coeffDict());
        CphiSscalar_.readIfPresent(coeffDict());
        CepsScalar_.readIfPresent(coeffDict());
        CpScalar_.readIfPresent(coeffDict());
        sigma_.readIfPresent(coeffDict());
        gN_.readIfPresent(coeffDict());
        maxK_.readIfPresent(coeffDict());
        ut_.readIfPresent(coeffDict());
        frictionalStressModel_->read();

        return true;
    }
    else
    {
        return false;
    }
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::SATFMdispersedModel::k() const
{
    dimensionedVector eSum
    (
        "eSum",
        dimensionSet(0, 0, 0, 0, 0, 0, 0),
        vector(1,1,1)
    );
    return k_&eSum;
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::SATFMdispersedModel::epsilon() const
{
    return Ceps_*pow(k(),3.0/2.0)/lm_;
}


Foam::tmp<Foam::volSymmTensorField>
Foam::RASModels::SATFMdispersedModel::R() const
{
    if (!anIsoTropicNut_) {
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
                2.0 * alpha_ * symm(R1_)
              - (nut_)*dev(twoSymm(fvc::grad(U_)))
            )
        );
    } else {
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
                2.0 * alpha_ * symm(R1_)
            )
        );
    }
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::SATFMdispersedModel::pPrime() const
{
    const volScalarField& rho = phase_.rho();
    tmp<volScalarField> tda(phase_.d());
    const volScalarField& da = tda();
    // Get strain rate tensor for frictional pressure models
    volTensorField gradU(fvc::grad(phase_.U()));
    boundGradU(gradU);
    volSymmTensorField D(symm(gradU));
    
    tmp<volScalarField> tpPrime
    (
        pos(alpha_ - alphaMinFriction_)
       *frictionalStressModel_->frictionalPressurePrime
        (
            phase_,
            alphaMinFriction_,
            alphaMax_,
            da,
            rho,
            dev(D)
        )
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
Foam::RASModels::SATFMdispersedModel::pPrimef() const
{
    return fvc::interpolate(pPrime());
}

Foam::tmp<Foam::volVectorField>
Foam::RASModels::SATFMdispersedModel::divStress() const
{
    const volScalarField& rho = phase_.rho();
    tmp<volScalarField> tda(phase_.d());
    const volScalarField& da = tda();
    // Get strain rate tensor for frictional pressure models
    volTensorField gradU(fvc::grad(phase_.U()));
    boundGradU(gradU);
    volSymmTensorField D(symm(gradU));
    
    volTensorField R1(R1_);
    boundStress(R1);
    R1.correctBoundaryConditions();
    
    tmp<volVectorField> tDivStress
    (
        pos(alpha_ - alphaMinFriction_)
       *fvc::grad
        (
            frictionalStressModel_->frictionalPressure
            (
                phase_,
                alphaMinFriction_,
                alphaMax_,
                da,
                rho,
                dev(D)
            )
         )
     + pos(alpha_ - residualAlpha_)
      *fvc::div
       (
           2.0
         * alpha_
         * rho_
         * R1
       )
    );

    return tDivStress;
}

Foam::tmp<Foam::volSymmTensorField>
Foam::RASModels::SATFMdispersedModel::devRhoReff() const
{
    if (!anIsoTropicNut_) {
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
                2.0 * alpha_ * rho_ * symm(R1_)
              - (rho_*nut_)*dev(twoSymm(fvc::grad(U_)))
            )
        );
    } else {
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
                2.0 * alpha_ * rho_ * symm(R1_)
            )
        );
    }
}


Foam::tmp<Foam::fvVectorMatrix>
Foam::RASModels::SATFMdispersedModel::divDevRhoReff
(
    volVectorField& U
) const
{
    if (!anIsoTropicNut_) {
        return
        pos(alpha_ - residualAlpha_)
      * (
          - fvm::laplacian(rho_*nut_, U)
          - fvc::div
            (
                (rho_*nut_)*dev2(T(fvc::grad(U)))
            )
         );
    } else {
        return
        pos(alpha_ - residualAlpha_)
      * (
          - fvm::laplacian(rho_*nuFric_, U)
          - fvc::div
            (
               (rho_*nuFric_)*dev2(T(fvc::grad(U)))
            )
         );
    }
}

void Foam::RASModels::SATFMdispersedModel::boundStress
(
    volTensorField& R
) const
{
    scalar RMin = -ut_.value()*ut_.value();
    scalar RMax = -RMin;

    R.max
    (
        dimensionedTensor
        (
            "zero",
            R.dimensions(),
            tensor
            (
                  0, 0.99*RMin, 0.99*RMin,
                  0.99*RMin, 0, 0.99*RMin,
                  0.99*RMin, 0.99*RMin, 0
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
                  RMax, 0.99*RMax, 0.99*RMax,
                  0.99*RMax, RMax, 0.99*RMax,
                  0.99*RMax, 0.99*RMax, RMax
            )
        )
    );
}

void Foam::RASModels::SATFMdispersedModel::boundNormalStress
(
    volVectorField& k
) const
{
    scalar kMin = 1.e-7;
    scalar kMax = maxK_.value();

    k.max
    (
        dimensionedVector
        (
            "zero",
            k.dimensions(),
            vector
            (
                kMin,
                kMin,
                kMin
            )
        )
     );

    k.min
    (
        dimensionedVector
        (
            "maxvalue",
            k.dimensions(),
            vector
            (
                kMax,
                kMax,
                kMax
            )
        )
    );
}

void Foam::RASModels::SATFMdispersedModel::boundxiPhiS
(
    volVectorField& xi
) const
{
    scalar xiMin = -0.99;
    scalar xiMax = 0.99;

    xi.max
    (
        dimensionedVector
        (
            "minXi",
            xi.dimensions(),
            vector
            (
                xiMin,
                xiMin,
                xiMin
            )
        )
    );
    xi.min
    (
        dimensionedVector
        (
            "maxXi",
            xi.dimensions(),
            vector
            (
                xiMax,
                xiMax,
                xiMax
            )
        )
    );
}

void Foam::RASModels::SATFMdispersedModel::boundCorrTensor
(
    volTensorField& R
) const
{
    scalar xiMin = -0.99;
    scalar xiMax = 0.99;

    R.max
    (
        dimensionedTensor
        (
            "zero",
            R.dimensions(),
            tensor
            (
                  xiMax, xiMin, xiMin,
                  xiMin, xiMax, xiMin,
                  xiMin, xiMin, xiMax
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
                  xiMax, xiMax, xiMax,
                  xiMax, xiMax, xiMax,
                  xiMax, xiMax, xiMax
            )
        )
    );
}

void Foam::RASModels::SATFMdispersedModel::boundGradU
(
    volTensorField& R
) const
{
    scalar sMin = -1.0e3;
    scalar sMax =  1.0e3;

    R.max
    (
        dimensionedTensor
        (
            "zero",
            R.dimensions(),
            tensor
            (
                  sMin, sMin, sMin,
                  sMin, sMin, sMin,
                  sMin, sMin, sMin
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
                  sMax, sMax, sMax,
                  sMax, sMax, sMax,
                  sMax, sMax, sMax
            )
        )
    );
}


void Foam::RASModels::SATFMdispersedModel::correct()
{
    // Local references
    const twoPhaseSystem& fluid = refCast<const twoPhaseSystem>(phase_.fluid());
    volScalarField alpha(min(max(alpha_, 1.e-7),alphaMax_));
    volScalarField alpha2 = 1.0 - alpha;
    const volScalarField& rho = phase_.rho();
    const surfaceScalarField& alphaRhoPhi = alphaRhoPhi_;
    const surfaceScalarField& phi1 = phi_;
    const volVectorField& U = U_;
    
    const volScalarField& rho2 = fluid.otherPhase(phase_).rho();
    
    // cont. Phase velocity
    const volVectorField& Uc_ = fluid.otherPhase(phase_).U();
    
    // slip velocity
    volVectorField uSlip = Uc_ - U;
    
    volVectorField UcZero
    (
        IOobject
        (
            "Uzero",
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Uc_,
        zeroGradientFvPatchField<vector>::typeName
     );
    UcZero.correctBoundaryConditions();
    
    tmp<volScalarField> tda(phase_.d());
    const volScalarField& da = tda();

    dimensionedScalar kSmall("kSmall", k_.dimensions(), 1.0e-6);
    dimensionedScalar uSmall("uSmall", U_.dimensions(), 1.0e-6);
    dimensionedScalar lSmall("lSmall", dimLength, 1.0e-6);
    
    dimensionedVector eX
    (
        "eX",
        dimensionSet(0, 0, 0, 0, 0, 0, 0),
        vector(1,0,0)
    );
    dimensionedVector eY
    (
        "eY",
        dimensionSet(0, 0, 0, 0, 0, 0, 0),
        vector(0,1,0)
    );
    dimensionedVector eZ
    (
        "eZ",
        dimensionSet(0, 0, 0, 0, 0, 0, 0),
        vector(0,0,1)
    );
    
    dimensionedVector eSum
    (
        "eSum",
        dimensionSet(0, 0, 0, 0, 0, 0, 0),
        vector(1,1,1)
    );
    
    dimensionedTensor zeroR
    (
        "eSum",
        dimensionSet(0, 2, -2, 0, 0, 0, 0),
        tensor(0,0,0,0,0,0,0,0,0)
    );
    
    volTensorField gradU(fvc::grad(U_));
    boundGradU(gradU);
    volSymmTensorField D(dev(symm(gradU)));
    // compute S_{ij}S_{ij} (no summation over j!!)
    volTensorField SijSij =  magSqr(gradU&eX)*(eX*eX)
                           + magSqr(gradU&eY)*(eY*eY)
                           + magSqr(gradU&eZ)*(eZ*eZ);
    volVectorField SijSijV =  ((SijSij&eX)&eSum)*eX
                            + ((SijSij&eY)&eSum)*eY
                            + ((SijSij&eZ)&eSum)*eZ;
    // gradient of solids volume fraction
    volVectorField gradAlpha  = fvc::grad(alpha);
    
    // get turbulent kinetic energy of continuous-phase
    const volVectorField& kC_(mesh_.lookupObject<volVectorField>
                                     ("k." + fluid.otherPhase(phase_).name()));
    
    // get alphaP2Mean
    const volScalarField& alphaP2Mean2_(mesh_.lookupObject<volScalarField>
                                        ("alphaP2Mean." + fluid.otherPhase(phase_).name()));
    volScalarField alphaP2MeanO = max(alphaP2Mean2_,alphaP2Mean_);

    // simple filter for local smoothing
    simpleFilter filterS(mesh_);
    
    // get drag coefficient
    volScalarField beta
    (
        fluid.lookupSubModel<dragModel>
        (
            phase_,
            fluid.otherPhase(phase_)
        ).Ki()
    );
    volScalarField betaA = beta/(rho);
    beta *= alpha;
    // compute total k
    volScalarField km(k_&eSum);
    km.max(kSmall.value());
    
    
    // local reference to deltaF
    volScalarField deltaF(deltaF_);
    // compute grid size
    const cellList& cells = mesh_.cells();
    forAll(cells,cellI)
    {
        scalar deltaMaxTmp = 0.0;
        const labelList& cFaces = mesh_.cells()[cellI];
        const point& centrevector = mesh_.cellCentres()[cellI];
        
        forAll(cFaces, cFaceI)
        {
            label faceI = cFaces[cFaceI];
            const point& facevector = mesh_.faceCentres()[faceI];
            scalar tmp = mag(facevector - centrevector);
            if (tmp > deltaMaxTmp)
            {
                deltaMaxTmp = tmp;
            }
        }
        deltaF[cellI] = 2*deltaMaxTmp;
    }
    
    volScalarField wD = wallDist(mesh_).y();
    
    // correction for cases w/o walls
    // (since wall distance is then negative)
    deltaF_ = neg(wD)*deltaF + pos(wD)*min(deltaF,wD);
    // compute mixing length
    lm_ =  neg(wD)*Cmu_*deltaF + pos(wD)*min(Cmu_*deltaF,CmuWScalar_*wD);
    // correction for cyclic patches
    {
        const fvPatchList& patches = mesh_.boundary();
        
        forAll(patches, patchi) {
            const fvPatch& curPatch = patches[patchi];

            if (isA<cyclicAMIFvPatch>(curPatch)||isA<cyclicACMIFvPatch>(curPatch)) {
                forAll(curPatch, facei) {
                    label celli = curPatch.faceCells()[facei];
                    deltaF_[celli] = deltaF[celli];
                    lm_[celli]     = Cmu_[celli]*deltaF[celli];
                }
            }
        }
    }
    deltaF_.max(lSmall.value());
    lm_.max(lSmall.value());
    
    // compute nut
    nut_ = alpha*sqrt(km)*lm_;
    if (dynamicAdjustment_) {
        volScalarField alphaf(filter_(alpha));
        alphaf.max(residualAlpha_.value());
        volScalarField alphafP2(filter_(sqr(alpha)));
        alphafP2.max(sqr(residualAlpha_.value()));
        volScalarField alphafP2Mean(alphafP2 - sqr(alphaf));
        alphafP2Mean.max(sqr(residualAlpha_.value()));
        volVectorField Uf(filter_(alpha*U)/alphaf);
        volScalarField alpha2f(1.0 - alphaf);
        volVectorField Ucf(filter_(alpha2*UcZero)/alpha2f);
        volScalarField aUU(filter_(alpha*magSqr(U))/alphaf - magSqr(Uf));
        
        //Compute xiPhiS
        volVectorField xiPhiSNom =   (
                                          filter_(alpha*U)
                                        - alphaf*filter_(U)
                                      );
        
        volScalarField tmpDenX = alphafP2Mean
                              * (
                                    filter_(alpha*sqr(U&eX)) / alphaf
                                  - sqr(Uf&eX)
                                );
        volScalarField tmpDenY = alphafP2Mean
                               * (
                                    filter_(alpha*sqr(U&eY)) / alphaf
                                  - sqr(Uf&eY)
                                 );
        volScalarField tmpDenZ = alphafP2Mean
                               * (
                                    filter_(alpha*sqr(U&eZ)) / alphaf
                                  - sqr(Uf&eZ)
                                );
       
        tmpDenX.max(SMALL);
        tmpDenY.max(SMALL);
        tmpDenZ.max(SMALL);

        xiPhiS_ =  eX
                 * (
                        filterS((xiPhiSNom&eX)*sqrt(tmpDenX))/filterS(tmpDenX)
                    )
                 + eY
                 * (
                        filterS((xiPhiSNom&eY)*sqrt(tmpDenY))/filterS(tmpDenY)
                    )
                 + eZ
                 * (
                        filterS((xiPhiSNom&eZ)*sqrt(tmpDenZ))/filterS(tmpDenZ)
                    );
        // limit xiPhiS_
        boundxiPhiS(xiPhiS_);
        
        // compute triple correlation
        volScalarField xiPhiGGnom
        (
            filter_(alpha*magSqr(UcZero))
          - alphaf*filter_(magSqr(UcZero))
          - 2.0*(Ucf&(filter_(alpha*UcZero) - alphaf*filter_(UcZero)))
        );
        volScalarField xiPhiGGden
        (
            sqrt(alphafP2Mean)
           *max(filter_(magSqr(UcZero)) - 2.0*(Ucf&filter_(UcZero)) + magSqr(Ucf),sqr(uSmall))
        );
        xiPhiGG_ = filterS(xiPhiGGnom*xiPhiGGden)/filterS(sqr(xiPhiGGden));
        // smooth and limit xiPhiGG_
        xiPhiGG_.max(-0.99);
        xiPhiGG_.min(0.99);
        
        // compute correlation coefficient between gas phase and solid phase velocity
        volScalarField xiGSnum
        (
            filter_(alpha*(UcZero&U))/alphaf
          - (filter_(alpha*UcZero) & Uf)/alphaf
        );
        volScalarField xiGSden
        (
            sqrt(max(filter_(alpha2*magSqr(UcZero))/alpha2f - magSqr(Ucf),kSmall))
           *sqrt(max(aUU,kSmall))
        );
        
        //xiGS_ = xiGSnum/xiGSden;
        xiGS_ = filterS(xiGSnum*xiGSden)/filterS(sqr(xiGSden));
                 
        // smooth and regularize xiGS_ (xiGS_ is positive)
        //xiGS_.max(0);
        xiGS_ = mag(xiGS_);
        xiGS_.min(sqrt(2.0));
        
        // xiPhiDivU
        volScalarField divU(fvc::div(U));
        volScalarField divUf(0.5*fvc::div(Uf));
        // volScalarField divUf(filter_(alpha*divU)/alphaf);
        divUf.max(SMALL);
        volScalarField divUfL(0.25*mag(fvc::laplacian(aUU)));
        divUfL.max(SMALL);
        volScalarField xiPhiDivUnum
        (
            filter_(alpha*divU)
          - alphaf*filter_(divU)
        );
        /*
        volScalarField xiPhiDivUden
        (
            sqrt(alphafP2Mean)
           *sqrt
            (
                max
                (
                    filter_(alpha*sqr(divU))/alphaf
                  - sqr(divUf)
                    ,
                    dimensionedScalar("small",dimensionSet(0,0,-2,0,0),1.0e-7)
                )
            )
        );
         */
        volScalarField xiPhiDivUden
        (
            sqrt(alphafP2Mean)
           *sqrt(divUfL)
        );
    
        xiPhiDivU_ = filterS(xiPhiDivUnum*xiPhiDivUden)/filterS(sqr(xiPhiDivUden));
        xiPhiDivU_.max(-1.0);
        xiPhiDivU_.min(1.0);
        
        //xiPhi2DivU
        volScalarField xiPhi2DivUnum
        (
            filter_(sqr(alpha)*divU)
          - alphaf
           *(
                2.0*filter_(alpha*divU)
              - alphaf*filter_(divU)
              - alphaf*divUf
            )
          - alphafP2*divUf
        );
        volScalarField xiPhi2DivUden
        (
            alphafP2Mean
           *sqrt(divUfL)
        );
        
        xiPhi2DivU_ = filterS(xiPhi2DivUnum*xiPhi2DivUden)/filterS(sqr(xiPhi2DivUden));
        xiPhi2DivU_.max(-1.0);
        xiPhi2DivU_.min(1.0);
        
        // compute mixing length dynamically
        /*
        volScalarField Lij  = filter_(alpha*magSqr(U))/alphaf - magSqr(Uf);
        Lij.max(SMALL);
        volScalarField magSqrDf = filter_(magSqr(D));
        magSqrDf.max(SMALL);
        volSymmTensorField Df   = filter_(D);
        volScalarField Mij = sqr(deltaF_)*(4.0*magSqr(Df) - magSqrDf);
        volScalarField MijMij = filterS(sqr(Mij));
        MijMij.max(SMALL);
        
        volScalarField CmuT = 0.5*(filterS(Lij * Mij)/(MijMij));
        
        CmuT.min(4.0*sqr(CmuScalar_.value()));
                
        Cmu_ = pos(alpha_ - residualAlpha_)*sqrt(0.5*(CmuT + mag(CmuT)) + scalar(1.0e-2))
             + neg(alpha_ - residualAlpha_)*CmuScalar_;
        */
        Cmu_    = CmuScalar_;
        
        // Currently no dynamic procedure for Ceps and Cp
        // Set Ceps
        // dynamic procedure for Ceps
        /*
        volScalarField LijEps = alpha*(nut_)*(magSqrDf - magSqr(Df));
        volScalarField MijEps = pow(alpha*Lij,1.5)/(2.0*deltaF_);
        volScalarField MijMijEps = filterS(sqr(MijEps));
        MijMijEps.max(SMALL);
        
        volScalarField CepsT = filterS(LijEps*MijEps)/(MijMijEps);
        
        Ceps_ = 0.5*pos(alpha_ - residualAlpha_)*(CepsT + mag(CepsT))
                  + neg(alpha_ - residualAlpha_);
        
        Ceps_.min(10.0);
        */
        Ceps_ = pos(alpha_ - residualAlpha_)*CepsScalar_
              + neg(alpha_- residualAlpha_);
        // compute CphiS
        CphiS_ = CphiSscalar_;
        // Set Cp
        Cp_     = CpScalar_;
    } else {
        volVectorField xiPhiSDir = uSlip
                                  /max(mag(uSlip),uSmall);
        xiPhiS_    = -(xiPhiSolidScalar_)*xiPhiSDir;
        xiPhiDivU_ = xiPhiDivUScalar_;
        xiPhiGG_   = scalar(0.0);
        xiGS_      = xiGSScalar_;
        Cmu_       = CmuScalar_;
        Ceps_      = CepsScalar_;
        Cp_        = CpScalar_;
        // compute CphiS
        CphiS_     = CphiSscalar_;
    }
    // xiPhiDivU_ = xiPhiDivUScalar_/alphaMax_*(scalar(1.0) - alpha/alphaMax_);
    // compute xiGatS
    xiGatS_ =  scalar(1.0) + xiPhiGG_*sqrt(alphaP2MeanO)
            / max(alpha*alpha2*(scalar(1.0) - xiPhiGG_*sqrt(alphaP2MeanO)/alpha2),residualAlpha_);
    xiGatS_.max(SMALL);
    xiGatS_.min(2.0);

    
    // Compute k_
    // ---------------------------
    if (!equilibriumK_) {
        volVectorField pDil = Cp_*alpha*(rho-rho2)*gN_*sqrt(2.0*alphaP2MeanO);
        
        volTensorField R1t(R1_);
        if (!anIsoTropicNut_) {
            R1t -= 0.5*nut_*dev(gradU + gradU.T());
        }
        // compute production term according to Reynolds-stress model
        volTensorField gradUR1 = 0.5*((R1t&gradU) + ((gradU.T())&(R1t.T())));
        
        // volTensorField gradUR1 = 0.5*((R1_&gradU) + (R1_.T()&gradU.T()));
        shearProd_ =  (
                         (gradUR1&&(eX*eX))*(eX)
                       + (gradUR1&&(eY*eY))*(eY)
                       + (gradUR1&&(eZ*eZ))*(eZ)
                      );
        // volScalarField coeffDissipation(Ceps_*alpha*rho/lm_);
        
        fv::options& fvOptions(fv::options::New(mesh_));
        
        // Construct the transport equation for k
        Info << "Solving k-equation (dispersed phase) ... " << endl;
        fvVectorMatrix kEqn
        (
            fvm::ddt(alpha, rho, k_)
          + fvm::div(alphaRhoPhi, k_)
          // + fvm::SuSp(-(fvc::ddt(alpha, rho) + fvc::div(alphaRhoPhi)), k_)
          // diffusion with anisotropic diffusivity
          - fvm::laplacian(alpha*rho*lm_
                                * (
                                     (sqrt(k_&eX)*(eX*eX))
                                   + (sqrt(k_&eY)*(eY*eY))
                                   + (sqrt(k_&eZ)*(eZ*eZ))
                                   )
                                 / (sigma_)
                           , k_
                           , "laplacian(kappa,k)"
                         )
          /*
          - fvm::laplacian(
                             alpha*rho*sqrt(km)*lm_/(sigma_),
                             k_,
                             "laplacian(kappa,k)"
                         )
          */
         ==
          // some source terms are explicit since fvm::Sp()
          // takes solely scalars as first argument.
          // ----------------
          // shear production
          - 2.0*alpha
               *rho
               *shearProd_
          // interfacial work (--> energy transfer)
          + 2.0*beta
               *(
                    xiGS_
                  * (
                        sqrt((kC_&eX)*(k_&eX))*eX
                      + sqrt((kC_&eY)*(k_&eY))*eY
                      + sqrt((kC_&eZ)*(k_&eZ))*eZ
                    )
                )
          + fvm::Sp(-2.0*beta,k_)
          // pressure dilation & dissipation
          // - (coeffDissipation*(k_&eX) + (pDil&eX)*(xiPhiS_&eX))*sqrt(k_&eX)*eX
          // - (coeffDissipation*(k_&eY) + (pDil&eY)*(xiPhiS_&eY))*sqrt(k_&eY)*eY
          // - (coeffDissipation*(k_&eZ) + (pDil&eZ)*(xiPhiS_&eZ))*sqrt(k_&eZ)*eZ
          - ((pDil&eX)*(xiPhiS_&eX))*sqrt(k_&eX)*eX
          - ((pDil&eY)*(xiPhiS_&eY))*sqrt(k_&eY)*eY
          - ((pDil&eZ)*(xiPhiS_&eZ))*sqrt(k_&eZ)*eZ
          // dissipation
          + fvm::Sp(-Ceps_*alpha*rho*sqrt(km)/deltaF_,k_)
          // + fvm::Sp(-Ceps_*alpha*rho*sqrt(D&&D),k_)
          + fvOptions(alpha, rho, k_)
        );

        kEqn.relax();
        fvOptions.constrain(kEqn);
        kEqn.solve();
        fvOptions.correct(k_);
    }
    else {
        // no dynamic adjustment for Ceps in case of equilibrium
        Ceps_   = CepsScalar_;
        // Equilibrium => dissipation == production
        // Schneiderbauer (2017), equ. (55)
        // pressure dilation is negligible in the solid phase
        forAll(cells,cellI)
        {
            for (int i=0; i<3; i++) {
                k_[cellI].component(i) =
                    sqr(
                         - betaA[cellI]*lm_[cellI]
                         + Foam::sqrt(
                              sqr(betaA[cellI]*lm_[cellI])
                            + 2.0 * lm_[cellI]
                            * Foam::max(
                                 lm_[cellI]*SijSijV[cellI].component(i)
                               + xiGS_[cellI]*betaA[cellI]*Foam::sqrt(Foam::max(kC_[cellI].component(i),kSmall.value()))
                                ,0.0
                              )
                           )
                        ) / sqr(Ceps_[cellI]);
            }
        }
        
        
    }
    // limit k_
    boundNormalStress(k_);
    // correct BCs
    k_.correctBoundaryConditions();

    //- compute variance of solids volume fraction
    // update km
    km = k_&eSum;
    km.max(kSmall.value());
    
    Info << "Computing nut (dispersed phase) ... " << endl;
    nut_ = alpha*sqrt(km)*lm_;
    
    // compute fields for transport equation for phiP2
    volScalarField divU(fvc::div(phi_));
    volScalarField dissPhiP2 = CphiS_ * Ceps_ * sqrt(km)/deltaF_;
    volScalarField denom = divU + dissPhiP2;
    volScalarField xiKgradAlpha = (
                                         ((sqrt(k_&eX) * (gradAlpha&eX) * (xiPhiS_&eX)))
                                       + ((sqrt(k_&eY) * (gradAlpha&eY) * (xiPhiS_&eY)))
                                       + ((sqrt(k_&eZ) * (gradAlpha&eZ) * (xiPhiS_&eZ)))
                                   );
    
    Info << "Computing alphaP2Mean (dispersed phase) ... " << endl;
    volScalarField alpha1(alpha);
    alpha1.min(0.99*alphaMax_.value());
    volScalarField alphaM(alphaMax_-alpha1);
    volScalarField g0(1.0/(1.0-alpha1/alphaMax_));
    
    volScalarField alphaL2
    (
        6.0*sqr(alpha)
       /(g0*(g0 + 2.0))
    );
    
    //volScalarField alphaL2(min(alpha1*alphaM,alphaL));
    alphaP2Mean_.max(SMALL);
    if (!equilibriumPhiP2_) {
        // Construct the transport equation for alphaP2Mean
        fvScalarMatrix phiP2Eqn
        (
            fvm::ddt(alphaP2Mean_)
          + fvm::div(phi1, alphaP2Mean_)
          // diffusion
          - fvc::div(
                        (
                           alpha*lm_
                         * (
                              (sqrt(k_&eX)*(eX*eX))
                            + (sqrt(k_&eY)*(eY*eY))
                            + (sqrt(k_&eZ)*(eZ*eZ))
                           )
                         / (sigma_)
                       )
                     & (fvc::grad(alphaP2Mean_/alpha))
                    )
          /*
          - fvc::div(
                       alpha*sqrt(km)*lm_/sigma_
                     * fvc::grad(alphaP2Mean_/alpha)
                    )
          
           */
           /*
          - fvm::laplacian(
                           sqrt(km)*lm_/(sigma_),
                           alphaP2Mean_
                         )
           */
         ==
          // production/dissipation
          - fvm::SuSp(divU,alphaP2Mean_)
          - fvm::SuSp(2.0*xiKgradAlpha/sqrt(alphaP2Mean_),alphaP2Mean_)
          - fvm::SuSp
            (
                (
                    2.0*xiPhiDivU_*alpha/sqrt(alphaP2Mean_)
                  - xiPhi2DivU_
                )
               *sqrt(mag(fvc::laplacian(k_)))
            ,
                alphaP2Mean_
            )
          //+ fvm::Sp(-dissPhiP2,alphaP2Mean_)
        );

        phiP2Eqn.relax();
        phiP2Eqn.solve();
    } else {
        denom.max(SMALL);
        alphaP2Mean_ =   8.0
                       * sqr(xiKgradAlpha)
                       / sqr(denom);
    }
    // limit alphaP2Mean
    alphaP2Mean_ = min(
                         alphaP2Mean_,
                         alphaL2
                      );
    alphaP2Mean_.max(VSMALL);
    alphaP2Mean_.correctBoundaryConditions();
    
    // use k_normal for nut in stress tensor
    volScalarField kT
    (
        1.5
       *(
            (k_&eSum)
          - mag(k_ & U_)
          / (mag(U_)+uSmall)
        )
    );
    kT.max(kSmall.value());
    nut_ = alpha*sqrt(kT)*lm_;
       
    if (anIsoTropicNut_) {
        volScalarField alphaf = filter_(alpha);
        alphaf.max(residualAlpha_.value());
        volVectorField Uf = filter_(alpha*U)/alphaf;
        
        // compute correlation coefficients
        volTensorField xiUUnom = (filter_(alpha*(U*U))/alphaf - Uf*Uf);
        volScalarField xiUUden = max(tr(xiUUnom),kSmall);
        
        xiUU_ = fvc::average(xiUUnom*xiUUden)/fvc::average(sqr(xiUUden));
        
        // limit correlation coefficients
        boundCorrTensor(xiUU_);
        
        xiUU_.correctBoundaryConditions();
        // compute Reynolds-stress tensor
        forAll(cells,cellI)
        {
            for (int i=0; i<3; i++) {
                for (int j=0; j<3; j++) {
                    if (i!=j) {
                        R1_[cellI].component(j+i*3) =  (xiUU_[cellI].component(j+i*3))
                                *sqrt(k_[cellI].component(i)*k_[cellI].component(j));
                    } else {
                        R1_[cellI].component(j+i*3) =  sqrt(k_[cellI].component(i)*k_[cellI].component(j));
                    }
                }
            }
        }
    } else {
        R1_  = (k_&eX)*(eX*eX) + (k_&eY)*(eY*eY) + (k_&eZ)*(eZ*eZ);
    }
    
    R1_.correctBoundaryConditions();
    
    // Frictional pressure
    pf_ = frictionalStressModel_->frictionalPressure
    (
        phase_,
        alphaMinFriction_,
        alphaMax_,
        da,
        rho,
        dev(D)
    );

    nuFric_ = frictionalStressModel_->nu
    (
        phase_,
        alphaMinFriction_,
        alphaMax_,
        pf_/rho,
        da,
        //strain-rate fluctuations --> Srivastrava (2003)
        dev(D)// + sqrt(km)/lm_*symmTensor::I
    );

    // BCs for nut_
    const fvPatchList& patches = mesh_.boundary();
    volScalarField::Boundary& nutBf = nut_.boundaryFieldRef();
    
    forAll(patches, patchi) {
        const fvPatch& curPatch = patches[patchi];

        if (isA<wallFvPatch>(curPatch)) {
            scalarField& nutw = nutBf[patchi];
            const vectorField& faceAreas
                = mesh_.Sf().boundaryField()[patchi];
            const scalarField& magFaceAreas
                = mesh_.magSf().boundaryField()[patchi];
            
            forAll(curPatch, facei) {
                label celli = curPatch.faceCells()[facei];
                nutw[facei] = lm_[celli]
                             *mag(k_[celli] - (k_[celli]&faceAreas[facei])*faceAreas[facei]/sqr(magFaceAreas[facei]));
            }
        }
    }

    // Limit viscosity and add frictional viscosity
    nut_.min(maxNut_);
    nuFric_ = min(nuFric_, maxNut_ - nut_);
    nuFric_.min(maxNut_);
    nut_ += nuFric_;
    
    volScalarField unity
    (
         IOobject
         (
              "unity",
              U.time().timeName(),
              U.mesh(),
              IOobject::NO_READ,
              IOobject::NO_WRITE
         ),
         U.mesh(),
         dimensionedScalar("unity", dimless, 1.0)
    );
    
    Info << "SA-TFM (dispersed Phase):" << nl
         << "    max(nut)         = " << max(nut_).value() << nl
         << "    max(nutFric)     = " << max(nuFric_).value() << nl
         << "    max(phiP2/phi2)  = " << max(alphaP2Mean_/sqr(alpha)).value() << nl
         << "    max(k1)          = " << max(k()).value() << nl
         << "    mean(k1x)        = " << fvc::domainIntegrate(alpha*(k_&eX)).value()
                                        /fvc::domainIntegrate(alpha).value()
                                      << nl
         << "    mean(k1y)        = " << fvc::domainIntegrate(alpha*(k_&eY)).value()
                                        /fvc::domainIntegrate(alpha).value()
                                      << nl
         << "    mean(k1z)        = " << fvc::domainIntegrate(alpha*(k_&eZ)).value()
                                        /fvc::domainIntegrate(alpha).value()
                                      << nl
         << "    mean(phiP2/phi2) = " << fvc::domainIntegrate(alphaP2Mean_/sqr(alpha)).value()
                                        /fvc::domainIntegrate(unity).value()
                                      << nl
         << endl;
}


// ************************************************************************* //

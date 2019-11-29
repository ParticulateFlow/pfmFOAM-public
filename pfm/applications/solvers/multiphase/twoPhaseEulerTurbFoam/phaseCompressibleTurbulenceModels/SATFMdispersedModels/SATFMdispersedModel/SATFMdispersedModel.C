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
        coeffDict_.lookupOrDefault<scalar>("maxNut",1000)
    ),

    xiPhiSolidScalar_
    (
        "xiPhiSolidScalar",
        dimensionSet(0,0,0,0,0),
        coeffDict_.lookupOrDefault<scalar>("xiPhiS",-0.1)
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
        coeffDict_.lookupOrDefault<scalar>("Cmu",0.25)
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
        zeroGradientFvPatchField<scalar>::typeName
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
            IOobject::NO_WRITE
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
            "deltaF",
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
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedTensor("zero", dimensionSet(0, 2, -2, 0, 0),
                           tensor(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0))
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
        zeroGradientFvPatchField<scalar>::typeName
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
        xiGSScalar_.readIfPresent(coeffDict());
        CmuScalar_.readIfPresent(coeffDict());
        CphiSscalar_.readIfPresent(coeffDict());
        CepsScalar_.readIfPresent(coeffDict());
        CpScalar_.readIfPresent(coeffDict());
        sigma_.readIfPresent(coeffDict());
        gN_.readIfPresent(coeffDict());
        maxK_.readIfPresent(coeffDict());
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
    NotImplemented;
    return nut_;
}


Foam::tmp<Foam::volSymmTensorField>
Foam::RASModels::SATFMdispersedModel::R() const
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
           2.0 * pos(alpha_ - residualAlpha_) * alpha_ * symm(R1_)
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::SATFMdispersedModel::pPrime() const
{
    const volScalarField& rho = phase_.rho();
    tmp<volScalarField> tda(phase_.d());
    const volScalarField& da = tda();
    // Get strain rate tensor for frictional pressure models
    tmp<volTensorField> tgradU(fvc::grad(phase_.U()));
    const volTensorField& gradU(tgradU());
    volSymmTensorField D(symm(gradU));

    tmp<volScalarField> tpPrime
    (
        frictionalStressModel_->frictionalPressurePrime
        (
            phase_,
            alphaMinFriction_,
            alphaMax_,
            da,
            rho,
            dev(D)
        )
      * pos(alpha_-alphaMinFriction_)
    );

    volScalarField::Boundary& bpPrime =
        tpPrime.ref().boundaryFieldRef();

    forAll(bpPrime, patchi)
    {
        if (!bpPrime[patchi].coupled())
        {
            bpPrime[patchi] == 0;
        }
    }

    return tpPrime;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::RASModels::SATFMdispersedModel::pPrimef() const
{
    return fvc::interpolate(pPrime());
}


Foam::tmp<Foam::volSymmTensorField>
Foam::RASModels::SATFMdispersedModel::devRhoReff() const
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
            2.0 * alpha_ * dev(symm(R1_))
        )
    );
}


Foam::tmp<Foam::fvVectorMatrix>
Foam::RASModels::SATFMdispersedModel::divDevRhoReff
(
    volVectorField& U
) const
{
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
    
    if (!anIsoTropicNut_) {
        return
        pos(alpha_ - residualAlpha_)*
        (
          - fvm::laplacian(rho_*nut_, U)
          - fvc::div
            (
                (rho_*nut_)*dev2(T(fvc::grad(U)))
            )
          + fvc::div
            (
                 2.0
               * alpha_
               * rho_
               * (
                     (R1_&&(eX*eX))*(eX*eX)
                   + (R1_&&(eY*eY))*(eY*eY)
                   + (R1_&&(eZ*eZ))*(eZ*eZ)
                 )
            )
        );
    } else {
        return
        pos(alpha_ - residualAlpha_)*
        (
          - fvm::laplacian(rho_*nut_, U)
          - fvc::div
            (
                (rho_*nut_)*dev2(T(fvc::grad(U)))
            )
          + fvc::div
            (
                 2.0
               * alpha_
               * rho_
               * (
                     R1_
                 )
            )
        );
    }
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
    scalar xiMin = -1.0;
    scalar xiMax = 1.0;

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
    volScalarField alpha(max(alpha_, 1.e-7));
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
    
    volTensorField gradU(fvc::grad(U_));
    boundGradU(gradU);
    volSymmTensorField D(dev(symm(gradU)));
    // compute S_{ij}S_{ij} (no summation over j!!)
    volTensorField SijSij =  magSqr(gradU&eX)*(eX*eX)
                           + magSqr(gradU&eY)*(eY*eY)
                           + magSqr(gradU&eZ)*(eZ*eZ);
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
    simpleFilterADM filterS(mesh_);
    
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
    volScalarField km  = (k_ & eSum);
    km.max(kSmall.value());
    
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
        deltaF_[cellI] = 2*deltaMaxTmp;
    }
    
    volScalarField wD = wallDist(mesh_).y();
    
    // correction for cases w/o walls
    // (since wall distance is then negative)
    deltaF_ = neg(wD)*deltaF_ + pos(wD)*min(deltaF_,2.0*wD);
    deltaF_.max(lSmall.value());
    
    if (dynamicAdjustment_) {
        volScalarField alphaf = filter_(alpha);
        alphaf.max(residualAlpha_.value());
        volScalarField alphafP2 = filter_(sqr(alpha));
        alphafP2.max(sqr(residualAlpha_.value()));
        volVectorField Uf = filter_(alpha*U)/alphaf;
        volScalarField alpha2f = 1.0 - alphaf;
        volVectorField Ucf = filter_(alpha2*Uc_)/alpha2f;
        volScalarField aUU = filter_(alpha*magSqr(U))/alphaf - magSqr(Uf);
        
        //Compute xiPhiS
        volVectorField xiPhiSNom =   (
                                          filter_(alpha*U)
                                        - alphaf*filter_(U)
                                      );
        volScalarField tmpA = alphafP2-sqr(alphaf);
        tmpA.max(ROOTVSMALL);
        volScalarField tmpK = filter_(alpha*magSqr(U)) / alphaf - magSqr(Uf);
        tmpK.max(ROOTVSMALL);
        /*
        volScalarField tmpDenX = tmpA
                              * (
                                    filter_(alpha*sqr(U&eX)) / alphaf
                                  - sqr(Uf&eX)
                                );
        volScalarField tmpDenY = tmpA
                               * (
                                    filter_(alpha*sqr(U&eY)) / alphaf
                                  - sqr(Uf&eY)
                                 );
        volScalarField tmpDenZ = tmpA
                               * (
                                    filter_(alpha*sqr(U&eZ)) / alphaf
                                  - sqr(Uf&eZ)
                                );
       
        tmpDenX.max(ROOTVSMALL);
        tmpDenY.max(ROOTVSMALL);
        tmpDenZ.max(ROOTVSMALL);

        xiPhiS_ =  eX
                 * (
                        filterS(sqrt(tmpDenX)*(xiPhiSNom&eX))/filterS(tmpDenX)
                    )
                 + eY
                 * (
                        filterS(sqrt(tmpDenY)*(xiPhiSNom&eY))/filterS(tmpDenY)
                    )
                 + eZ
                 * (
                        filterS(sqrt(tmpDenZ)*(xiPhiSNom&eZ))/filterS(tmpDenZ)
                    );
        */
        xiPhiS_ = sqrt(3.0)*xiPhiSNom/sqrt(tmpA*tmpK);
        // smooth xiPhiS_
        xiPhiS_ = filterS(xiPhiS_);
        // limit xiPhiS_
        boundxiPhiS(xiPhiS_);
        
        // compute triple correlation
        volScalarField xiPhiGGnom =  filter_(alpha*magSqr(Uc_))
                                   - alphaf*filter_(magSqr(Uc_))
                                   - 2.0*(Ucf&(filter_(alpha*Uc_) - alphaf*filter_(Uc_)));
        volScalarField xiPhiGGden =  sqrt(max(alphafP2-sqr(alphaf),sqr(residualAlpha_)))
                                   * max(filter_(magSqr(Uc_))- 2.0*(Ucf&filter_(Uc_)) + magSqr(Ucf),sqr(uSmall));
        xiPhiGG_ = xiPhiGGnom/xiPhiGGden;
        // smooth and limit xiPhiGG_
        xiPhiGG_.max(-0.99);
        xiPhiGG_.min(0.99);
        xiPhiGG_ = filterS(xiPhiGG_);
        
        // compute correlation coefficient between gas phase and solid phase velocity
        volScalarField xiGSnum  = filter_(alpha*(Uc_&U))/alphaf - (filter_(alpha*Uc_) & Uf)/alphaf;
        volScalarField xiGSden  = sqrt(max(filter_(alpha*magSqr(Uc_))/alphaf-2.0*((filter_(alpha*Uc_)/alphaf)&Ucf)+magSqr(Ucf),kSmall))
                               * sqrt(max(aUU,kSmall));
        /*
        volScalarField xiGSnumP(xiGSnum*xiGSden);
        
        //set xiGSnumP to 0 at boundaries
        const fvPatchList& patches = mesh_.boundary();

        volScalarField::Boundary& xiGSnumPBf = xiGSnumP.boundaryFieldRef();

        forAll(patches, patchi)
        {
            if (!patches[patchi].coupled())
            {
                xiGSnumPBf[patchi] *= 0.0;
            }
        }
        */
        xiGS_ = xiGSnum/xiGSden;

        // smooth and regularize xiGS_ (xiGS_ is positive)
        xiGS_.max(-0.99);
        xiGS_.min(0.99);
        xiGS_ = filterS(xiGS_);
    
        // compute mixing length dynamically
        /*
        volScalarField Lij  = filter_(alpha*magSqr(U))/alphaf - magSqr(Uf);
        Lij.max(0);
        volScalarField Mij = sqr(deltaF_)*(4.0*magSqr(filter_(alpha*D)/alphaf) - filter_(alpha*magSqr(D))/alphaf);
        volScalarField MijMij = filterS(sqr(Mij));
        MijMij.max(SMALL);
        
        volScalarField CmuT = 0.5*(filterS(Lij * Mij)/(MijMij));
        
        CmuT.min(2.0*sqr(CmuScalar_.value()));
        CmuT.max(0.01);
        
        Cmu_ = pos(alpha_ - residualAlpha_)*sqrt(CmuT)
             + neg(alpha_- residualAlpha_)*CmuScalar_;
        */
        Cmu_    = CmuScalar_;
        
        // Currently no dynamic procedure for Ceps and Cp
        // Set Ceps
        Ceps_ = pos(alpha_ - residualAlpha_)*CepsScalar_
              + neg(alpha_- residualAlpha_);
        // compute CphiS
        CphiS_ = CphiSscalar_*Cmu_;
        // Set Cp
        Cp_     = CpScalar_;
    } else {
        volVectorField xiPhiSDir = uSlip
                                  /max(mag(uSlip),uSmall);
        xiPhiS_ = -(xiPhiSolidScalar_)*xiPhiSDir;
        xiPhiGG_ = scalar(0.0);
        xiGS_   = xiGSScalar_;
        Cmu_    = CmuScalar_;
        Ceps_   = CepsScalar_;
        Cp_     = CpScalar_;
        // compute CphiS
        CphiS_ = CphiSscalar_*CmuScalar_;
    }
    
    // compute mixing length
    lm_ = Cmu_*deltaF_;

    // compute xiGatS
    xiGatS_ =  scalar(1.0) + xiPhiGG_*sqrt(alphaP2MeanO)
            / max(alpha*alpha2*(scalar(1.0) - xiPhiGG_*sqrt(alphaP2MeanO)/alpha2),residualAlpha_);
    xiGatS_.max(1.0e-7);
    xiGatS_.min(2.0);
    // correct xiGS_
    xiGS_ *= sqrt(xiGatS_);
    
    // Compute k_
    // ---------------------------
    if (!equilibriumK_) {
        volVectorField pDil = Cp_*alpha*(rho-rho2)*gN_*sqrt(2.0*alphaP2MeanO);
        
        /*
        volTensorField R1t(R1_);
        if (!anIsoTropicNut_) {
            R1t -= nut_*dev(gradU + gradU.T());
        }
        // compute production term according to Reynolds-stres model
        volTensorField gradUR1 = 0.5*((R1t&gradU) + (R1t.T()&gradU.T()));
        */
        volTensorField gradUR1 = 0.5*((R1_&gradU) + (R1_.T()&gradU.T()));
        shearProd_ =   (gradUR1&&(eX*eX))*(eX)
                     + (gradUR1&&(eY*eY))*(eY)
                     + (gradUR1&&(eZ*eZ))*(eZ);
        
        fv::options& fvOptions(fv::options::New(mesh_));
        
        // Construct the transport equation for k
        Info << "Solving k-equation (dispersed phase) ... " << endl;
        fvVectorMatrix kEqn
        (
            fvm::ddt(alpha, rho, k_)
          + fvm::div(alphaRhoPhi, k_)
          - fvc::Sp((fvc::ddt(alpha, rho) + fvc::div(alphaRhoPhi)), k_)
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
          //- (coeffDissipation*(k_&eX) + (pDil&eX)*(xiPhiS_&eX))*sqrt(k_&eX)*eX
          //- (coeffDissipation*(k_&eY) + (pDil&eY)*(xiPhiS_&eY))*sqrt(k_&eY)*eY
          //- (coeffDissipation*(k_&eZ) + (pDil&eZ)*(xiPhiS_&eZ))*sqrt(k_&eZ)*eZ
          - ((pDil&eX)*(xiPhiS_&eX))*sqrt(k_&eX)*eX
          - ((pDil&eY)*(xiPhiS_&eY))*sqrt(k_&eY)*eY
          - ((pDil&eZ)*(xiPhiS_&eZ))*sqrt(k_&eZ)*eZ
          // dissipation
          + fvm::Sp(-Ceps_*alpha*rho*sqrt(km)/lm_,k_)
          // + fvm::Sp(-Ceps_*alpha*rho*sqrt(D&&D),k_)
          + fvOptions(alpha, rho, k_)
        );

        kEqn.relax();
        fvOptions.constrain(kEqn);
        kEqn.solve();
        fvOptions.correct(k_);
    }
    else {
        volVectorField SijSijV =  ((SijSij&eX)&eSum)*eX
                                + ((SijSij&eY)&eSum)*eY
                                + ((SijSij&eZ)&eSum)*eZ;
                                  
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
    
    // limit k before computing Reynolds-stresses
    boundNormalStress(k_);
    k_.correctBoundaryConditions();

    //- compute variance of solids volume fraction
    // update km
    km  = (k_ & eSum);
    km.max(kSmall.value());
    
    // compute fields for transport equation for phiP2
    volScalarField divU(fvc::div(U));
    volScalarField dissPhiP2 = CphiS_ * Ceps_ * sqrt(km)/lm_;
    volScalarField denom = mag(divU) + dissPhiP2;
    denom.max(SMALL);
    volScalarField xiKgradAlpha = (
                                 ((sqrt(k_&eX) * (gradAlpha&eX) * (xiPhiS_&eX)))
                               + ((sqrt(k_&eY) * (gradAlpha&eY) * (xiPhiS_&eY)))
                               + ((sqrt(k_&eZ) * (gradAlpha&eZ) * (xiPhiS_&eZ)))
                               );
    
    Info << "Computing alphaP2Mean (dispersed phase) ... " << endl;
    if (!equilibriumPhiP2_) {
        // Construct the transport equation for alphaP2Mean
        fvScalarMatrix phiP2Eqn
        (
            fvm::ddt(alphaP2Mean_)
          + fvm::div(phi1, alphaP2Mean_)
          - fvm::laplacian(
                           lm_
                         * (
                              (sqrt(k_&eX)*(eX*eX))
                            + (sqrt(k_&eY)*(eY*eY))
                            + (sqrt(k_&eZ)*(eZ*eZ))
                           )
                         / (sigma_)
                       ,
                         alphaP2Mean_
                       )
         ==
          // some source terms are explicit since fvm::Sp()
          // takes solely scalars as first argument.
          // ----------------
          // shear production
          - fvm::SuSp(
                         2.0*xiKgradAlpha
                       / sqrt(alphaP2Mean_+dimensionedScalar("small",dimensionSet(0,0,0,0,0), 1.0e-8))
                     ,alphaP2Mean_)
          - fvm::SuSp(divU,alphaP2Mean_)
          // production/dissipation
          + fvm::Sp(-dissPhiP2,alphaP2Mean_)
        );

        phiP2Eqn.relax();
        phiP2Eqn.solve();
    } else {
        alphaP2Mean_ =   8.0
                       * sqr(xiKgradAlpha)
                       / sqr(denom);
    }
    // limit alphaP2Mean_
    volScalarField alphaM = alphaMax_ - alpha;
    alphaM.max(0.0);
    volScalarField alphaL2 = sqr(min(alpha,alphaM));
    alphaP2Mean_ = min(
                         alphaP2Mean_,
                         0.99*alphaL2
                      );
    alphaP2Mean_.max(ROOTVSMALL);
    alphaP2Mean_.correctBoundaryConditions();
    
    {
        volScalarField alphaf = filter_(alpha);
        alphaf.max(residualAlpha_.value());
        volVectorField Uf = filter_(alpha*U)/alphaf;
        
        // compute correlation coefficients
        volTensorField xiUUnom = filter_(alpha*(U*U))/alphaf - Uf*Uf;
        volVectorField xiUUden = (
                                    sqrt(max(filter_(alpha*magSqr(U&eX))/alphaf - magSqr(Uf&eX),kSmall))*eX
                                  + sqrt(max(filter_(alpha*magSqr(U&eY))/alphaf - magSqr(Uf&eY),kSmall))*eY
                                  + sqrt(max(filter_(alpha*magSqr(U&eZ))/alphaf - magSqr(Uf&eZ),kSmall))*eZ
                                 );

        forAll(cells,cellI)
        {
            for (int i=0; i<3; i++) {
                for (int j=0; j<3; j++) {
                    xiUU_[cellI].component(j+i*3) = xiUUnom[cellI].component(j+i*3)
                                                  / (xiUUden[cellI].component(i)*xiUUden[cellI].component(j));
                }
            }
        }
        // limit correlation coefficients
        boundCorrTensor(xiUU_);
        xiUU_ = filterS(xiUU_);
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
        nut_ = 0.5*alpha*sqrt(dev(R1_)&&dev(R1_))/(sqrt(D&&D)+dimensionedScalar("small",dimensionSet(0,0,-1,0,0),SMALL));
        
        //set R1 to 0 at boundaries
        const fvPatchList& patches = mesh_.boundary();
        volTensorField::Boundary& R1Bf = R1_.boundaryFieldRef();

        forAll(patches, patchi)
        {
            if (patches[patchi].type() == "wall")
            {
                for (int i=0; i<3; i++) {
                    for (int j=0; j<3; j++) {
                        R1Bf[patchi].component(j+i*3) = 0;
                    }
                }
            }
        }
    }
    // R1_.correctBoundaryConditions();
  
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
        dev(D)
    );
    
    // Limit viscosity and add frictional viscosity
    nut_.min(maxNut_);
    nut_.correctBoundaryConditions();
    nuFric_ = min(nuFric_, maxNut_ - nut_);
    nuFric_.min(maxNut_);
    nut_ += nuFric_;
    
    Info << "SA-TFM (dispersed Phase):" << nl
         << "    max(nut)        = " << max(nut_).value() << nl
         << "    max(nutFric)    = " << max(nuFric_).value() << nl
         << "    max(phiP2/phi2) = " << max(alphaP2Mean_/sqr(alpha)).value() << nl
         << "    max(k1)         = " << max(k_&eSum).value() << nl
         << "    mean(k1x)       = " << fvc::domainIntegrate(alpha*(k_&eX)).value()
                                        /fvc::domainIntegrate(alpha).value()
                                     << nl
         << "    mean(k1y)       = " << fvc::domainIntegrate(alpha*(k_&eY)).value()
                                        /fvc::domainIntegrate(alpha).value()
                                     << nl
         << "    mean(k1z)       = " << fvc::domainIntegrate(alpha*(k_&eZ)).value()
                                        /fvc::domainIntegrate(alpha).value()
                                     << nl
         << "    mean(phi1P2)    = " << fvc::domainIntegrate(alpha*(alphaP2Mean_)).value()
                                        /fvc::domainIntegrate(alpha).value()
         << endl;
}


// ************************************************************************* //

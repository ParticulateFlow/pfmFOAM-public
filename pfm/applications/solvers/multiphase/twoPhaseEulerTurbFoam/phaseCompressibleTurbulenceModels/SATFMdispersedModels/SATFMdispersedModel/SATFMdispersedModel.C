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

    equilibrium_(coeffDict_.lookup("equilibrium")),
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
        // Set Boundary condition
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
        // Set Boundary condition
        zeroGradientFvPatchField<scalar>::typeName
    ),

    alphaP2Mean_
    (
        IOobject
        (
            IOobject::groupName("alphaP2Mean", phase.name()),
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 1.0e-7),
        // Set Boundary condition
        zeroGradientFvPatchField<scalar>::typeName
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
        dimensionedScalar("value", dimensionSet(0, 0, 0, 0, 0), 1.0e-2),
        // Set Boundary condition
        fixedValueFvPatchField<scalar>::typeName
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
        dimensionedScalar("value", dimensionSet(0, 0, 0, 0, 0), 0.4),
        // Set Boundary condition
        zeroGradientFvPatchField<scalar>::typeName
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
        dimensionedScalar("value", dimensionSet(0, 0, 0, 0, 0), 0.3),
        // Set Boundary condition
        zeroGradientFvPatchField<scalar>::typeName
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
                           tensor(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)),
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
        coeffDict().lookup("equilibrium") >> equilibrium_;
        coeffDict().lookup("dynamicAdjustment") >> dynamicAdjustment_;
        coeffDict().lookup("anIsoTropicNut") >> anIsoTropicNut_;
        alphaMax_.readIfPresent(coeffDict());
        alphaMinFriction_.readIfPresent(coeffDict());
        xiPhiSolidScalar_.readIfPresent(coeffDict());
        CmuScalar_.readIfPresent(coeffDict());
        CphiSscalar_.readIfPresent(coeffDict());
        CepsScalar_.readIfPresent(coeffDict());
        CpScalar_.readIfPresent(coeffDict());
        sigma_.readIfPresent(coeffDict());
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
          - (nut_)*dev(twoSymm(fvc::grad(U_)))
          + 2.0 * pos(alpha_ - residualAlpha_) * alpha_ *
            symm(R1_)
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
        )*pos(alpha_-alphaMinFriction_)
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
          - (rho_*nut_)*dev(twoSymm(fvc::grad(U_)))
          + 2.0 * pos(alpha_ - residualAlpha_) * alpha_ * rho_ *
            symm(R1_)
        )
    );
}


Foam::tmp<Foam::fvVectorMatrix>
Foam::RASModels::SATFMdispersedModel::divDevRhoReff
(
    volVectorField& U
) const
{
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
             2.0 * alpha_ * rho_ * R1_
        )
       *pos(alpha_ - residualAlpha_)
     
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
    scalar kMin = -1.0;
    scalar kMax = 1.0;

    R.max
    (
        dimensionedTensor
        (
            "zero",
            R.dimensions(),
            tensor
            (
                  kMax, kMin, kMin,
                  kMin, kMax, kMin,
                  kMin, kMin, kMax
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

void Foam::RASModels::SATFMdispersedModel::boundS
(
    volTensorField& R
) const
{
    scalar sMin = 1.0e-7;
    scalar sMax = 100.;

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
    volScalarField alpha(max(alpha_, residualAlpha_));
    const volScalarField& rho = phase_.rho();
    const surfaceScalarField& alphaRhoPhi = alphaRhoPhi_;
    const volVectorField& U = U_;
    
    // const volScalarField& rho2 = fluid.otherPhase(phase_).rho();
    
    // cont. Phase velocity
    const volVectorField& Uc_ = fluid.otherPhase(phase_).U();
    
    // slip velocity
    volVectorField uSlip = Uc_ - U;
    
    // gravity vector
    // const uniformDimensionedVectorField& g = mesh_.lookupObject<uniformDimensionedVectorField>("g");
    
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
    
    tmp<volTensorField> tgradU(fvc::grad(U_));
    const volTensorField& gradU(tgradU());
    volSymmTensorField D(dev(symm(gradU)));
    // compute S_{ij}S_{ij} (no summation over j!!)
    volTensorField SijSij =  magSqr(gradU&eX)*(eX*eX)
                           + magSqr(gradU&eY)*(eY*eY)
                           + magSqr(gradU&eZ)*(eZ*eZ);
    boundS(SijSij);
    // Set SijSij to 0 in cell with low alpha
    SijSij *= pos(alpha - residualAlpha_);
    
    // gradient of solids volume fraction
    volVectorField gradAlpha  = fvc::grad(alpha);
    
    // get turbulent kinetic energy of continuous-phase
    const volVectorField& kC_(mesh_.lookupObject<volVectorField>
                                     ("k." + fluid.otherPhase(phase_).name()));
    
    // get correlation coefficient between continuous and solid phase velocities
    const volScalarField& xiGS_(mesh_.lookupObject<volScalarField>
                             ("xiGS"));
    
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
    volScalarField km  = k_ & eSum;
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
        // compute xiPhiS
        volVectorField Uf = filter_(alpha*U)/alphaf;
        volScalarField aUU = filter_(alpha*(U&U))/alphaf - magSqr(Uf);


        volVectorField xiPhiSNom =   (
                                          filter_(alpha*U)
                                        - alphaf*filter_(U)
                                      );

        volScalarField xiPhiSDenomSqr =   (filter_(sqr(alpha))-sqr(alphaf))
                                        * (filter_(alpha*magSqr(U))/alphaf - magSqr(Uf));
        xiPhiSDenomSqr.max(kSmall.value());
        xiPhiS_ = sqrt(3.0)*xiPhiSNom/sqrt(xiPhiSDenomSqr);
       


        // smooth correlation coefficient
        xiPhiS_ = filterS(xiPhiS_);
        xiPhiS_ = 0.5*(
                        - mag(xiPhiS_)*gradAlpha
                              /(mag(gradAlpha)+dimensionedScalar("small",dimensionSet(0,-1,0,0,0),1.e-7))
                        + xiPhiS_
                      );
        //  xiPhiS_ = 0.5*(-mag(xiPhiS_)*uSlip/(mag(uSlip)+uSmall) + xiPhiS_);
        boundxiPhiS(xiPhiS_);
        
        // Currently no dynamic procedure for Ceps and Cp
        // Set Ceps
        Ceps_   = CepsScalar_;
        // compute CphiS
        CphiS_ = CphiSscalar_/CepsScalar_;
        // Set Cp
        Cp_     = CpScalar_;
        
        // compute mixing length dynamically
        volScalarField Lij  = filter_(alpha*magSqr(U))/alphaf - magSqr(Uf);
        Lij.max(0);
        volScalarField alpha2 = 1.0 - alpha;
        volScalarField alpha2f = 1.0 - alphaf;
        volVectorField Ucf = filter_(alpha2*Uc_)/alpha2f;
        volScalarField Lijc = filter_(alpha2*magSqr(Uc_))/alpha2f - magSqr(Ucf);
        Lijc.max(0);
        //Lij -= max(xiGS_*(sqrt(Lij*Lijc) - Lij),kSmall);
        //Lij.max(0);
        // volScalarField Mij = sqr(deltaF_)*(2.0*magSqr(dev(symm(fvc::grad(Uf)))) - filter_(magSqr(D)));
        volScalarField Mij = sqr(deltaF_)*(4.0*magSqr(filter_(alpha*D)/alphaf) - filter_(alpha*magSqr(D))/alphaf);
        volScalarField MijMij = filterS(Mij * Mij);
        MijMij.max(SMALL);
        volScalarField CmuT = 0.5*filterS(Lij * Mij)/(MijMij);
        
        CmuT = 0.5*(mag(CmuT) + CmuT);
        CmuT.max(0);
        
        Cmu_ = sqrt(CmuT);
        Cmu_ = filterS(Cmu_);

        Cmu_.min(2.0*CmuScalar_.value());
        Cmu_.max(0.01*CmuScalar_.value());
        
    } else {
        xiPhiS_ = - xiPhiSolidScalar_*gradAlpha
                   /max(mag(gradAlpha),dimensionedScalar("small",dimensionSet(0,-1,0,0,0),1.e-7));
        Cmu_    = CmuScalar_;
        Ceps_   = CepsScalar_;
        Cp_     = CpScalar_;
        // compute CphiS
        CphiS_ = CphiSscalar_/CepsScalar_;
    }
    
    // compute mixing length
    lm_ = Cmu_*deltaF_;

    
    // Compute k_
    // ---------------------------
    if (!equilibrium_) {
        // volVectorField pDil = Cp_*alpha*(rho-rho2)*g*sqrt(2.0*alphaP2MeanO);
        
        fv::options& fvOptions(fv::options::New(mesh_));
        
        // Construct the transport equation for k
        // --> Stefanie
        Info << "Solving k-equation (dispersed phase) ... " << endl;
        fvVectorMatrix kEqn
        (
            fvm::ddt(alpha, rho, k_)
          + fvm::div(alphaRhoPhi, k_)
          - fvc::Sp(fvc::ddt(alpha, rho) + fvc::div(alphaRhoPhi), k_)
          // diffusion with anisotropic diffusivity
          - fvm::laplacian(alpha*rho*lm_
                                * (
                                     (sqrt(k_&eX)*(eX*eX))
                                   + (sqrt(k_&eY)*(eY*eY))
                                   + (sqrt(k_&eZ)*(eZ*eZ))
                                   )
                                 / (sigma_)
                           , k_, "laplacian(kappa,k)")
         ==
          // some source terms are explicit since fvm::Sp()
          // takes solely scalars as first argument.
          // ----------------
          // shear production
            2.0*lm_
               *alpha
               *rho
               *(
                    (((SijSij&eX)&eSum)*(k_&eX))*eX
                  + (((SijSij&eY)&eSum)*(k_&eY))*eY
                  + (((SijSij&eZ)&eSum)*(k_&eZ))*eZ
                )
               /sqrt(km)
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
          // pressure dilation
          // - ((pDil&eX)*(xiPhiS_&eX)*sqrt(k_&eX))*eX
          // - ((pDil&eZ)*(xiPhiS_&eY)*sqrt(k_&eY))*eY
          // - ((pDil&eY)*(xiPhiS_&eZ)*sqrt(k_&eZ))*eZ
          // dissipation
          + fvm::Sp(-Ceps_*alpha*rho*sqrt(km)/lm_,k_)
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
    km  = k_ & eSum;
    km.max(kSmall.value());
    volScalarField divU(fvc::div(U));
    volScalarField denom = divU + CphiS_ * Ceps_ * sqrt(km)/lm_;
    denom.max(kSmall.value());
    
    Info << "Computing alphaP2Mean (dispersed phase) ... " << endl;
    if (dynamicAdjustment_) {
        volScalarField xiKgradAlpha = (
                                         ((sqrt(k_&eX) * (gradAlpha&eX) * (xiPhiS_&eX)))
                                       + ((sqrt(k_&eY) * (gradAlpha&eY) * (xiPhiS_&eY)))
                                       + ((sqrt(k_&eZ) * (gradAlpha&eZ) * (xiPhiS_&eZ)))
                                       );
        alphaP2Mean_ =   8.0
                       * sqr(xiKgradAlpha)*neg(xiKgradAlpha)
                       / sqr(denom);
    } else {
        alphaP2Mean_ =   8.0
                       * magSqr(xiPhiS_)
                       * sqr(
                                (sqrt(k_&eX) * mag(gradAlpha&eX))
                              + (sqrt(k_&eY) * mag(gradAlpha&eY))
                              + (sqrt(k_&eZ) * mag(gradAlpha&eZ))
                         )
                       / sqr(denom);
    }
    // limti alphaP2Mean_
    alphaP2Mean_.max(sqr(residualAlpha_.value()));
    volScalarField alphaM = alphaMax_ - alpha;
    alphaM.max(0.0);
    alphaP2Mean_ = min(
                         alphaP2Mean_,
                         alpha*alphaM/(alphaMax_*alphaMax_)
                      );

    // compute nut_ (Schneiderbauer, 2017; equ. (34))
    nut_ = pos(alpha - residualAlpha_)*alpha*sqrt(km)*lm_;
    if (anIsoTropicNut_) {
        nut_ *= 0.;
        volScalarField alphaf = filter_(alpha);
        alphaf.max(residualAlpha_.value());
        volVectorField Uf = filter_(alpha*U)/alphaf;
        // compute correlation coefficients
        xiUU_ = 3.0*(filter_(alpha*(U*U))/alphaf - Uf*Uf)
                /max(filter_(alpha*(U&U))/alphaf - (Uf&Uf),kSmall);
        // limit and smooth correlation coefficients
        xiUU_ = filterS(xiUU_);
        boundCorrTensor(xiUU_);
        
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
        R1_ = (k_&eX)*(eX*eX) + (k_&eY)*(eY*eY) + (k_&eZ)*(eZ*eZ);
    }
    
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
    nuFric_ = min(nuFric_, maxNut_ - nut_);
    nuFric_.min(maxNut_);
    nut_ += nuFric_;
    
    Info<< "SA-TFM (dispersed Phase):" << nl
        << "    max(nut) = " << max(nut_).value() << nl
        << "    max(nutFric) = " << max(nuFric_).value() << nl
        << "    max(k_) = " << max(k_&eSum).value() << endl;
}


// ************************************************************************* //

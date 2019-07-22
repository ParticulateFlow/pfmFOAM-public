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

#include "SATFMcontinuousModel.H"
#include "mathematicalConstants.H"
#include "twoPhaseSystem.H"
#include "wallDist.H"
#include "simpleFilter.H"
#include "uniformDimensionedFields.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RASModels::SATFMcontinuousModel::SATFMcontinuousModel
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

    equilibrium_(coeffDict_.lookup("equilibrium")),
    dynamicAdjustment_(coeffDict_.lookup("dynamicAdjustment")),

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

    xiPhiContScalar_
    (
        "xiPhiContScalar",
        dimensionSet(0,0,0,0,0),
        coeffDict_.lookupOrDefault<scalar>("xiPhiG",-0.5)
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
        coeffDict_.lookupOrDefault<scalar>("Cmu",0.4)
    ),

    CphiGscalar_
    (
        "CphiGscalar",
        dimensionSet(0,0,0,0,0),
        coeffDict_.lookupOrDefault<scalar>("CphiG",0.4)
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
        "sigmaC",
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

    xiPhiG_
    (
        IOobject
        (
            "xiPhiG",
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
        dimensionedVector("value", dimensionSet(0, 0, 0, 0, 0), vector(-0.5,-0.5,-0.5)),
        // Set Boundary condition
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
            IOobject::NO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("value", dimensionSet(0, 0, 0, 0, 0), 1.0),
        // Set Boundary condition
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
            IOobject::NO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("value", dimensionSet(0, 0, 0, 0, 0), 0.9),
        // Set Boundary condition
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
            IOobject::NO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("value", dimensionSet(0, 0, 0, 0, 0), 1.0),
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
            IOobject::NO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("value", dimensionSet(0, 0, 0, 0, 0), 0.4),
        // Set Boundary condition
        zeroGradientFvPatchField<scalar>::typeName
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
        dimensionedScalar("value", dimensionSet(0, 0, 0, 0, 0), 1.0),
        // Set Boundary condition
        zeroGradientFvPatchField<scalar>::typeName
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
            "lm",
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("value", dimensionSet(0, 1, 0, 0, 0), 1.e-2)
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

Foam::RASModels::SATFMcontinuousModel::~SATFMcontinuousModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::RASModels::SATFMcontinuousModel::read()
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

        xiPhiContScalar_.readIfPresent(coeffDict());
        xiGSScalar_.readIfPresent(coeffDict());
        CmuScalar_.readIfPresent(coeffDict());
        CphiGscalar_.readIfPresent(coeffDict());
        CepsScalar_.readIfPresent(coeffDict());
        CpScalar_.readIfPresent(coeffDict());
        sigma_.readIfPresent(coeffDict());
        maxK_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::SATFMcontinuousModel::k() const
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
Foam::RASModels::SATFMcontinuousModel::epsilon() const
{
    NotImplemented;
    return nut_;
}


Foam::tmp<Foam::volSymmTensorField>
Foam::RASModels::SATFMcontinuousModel::R() const
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
          + 2.0 * alpha_ * rho_ *
            symm(
                     (k_&eX)*(eX*eX)
                   + (k_&eY)*(eY*eY)
                   + (k_&eZ)*(eZ*eZ)
                 )
        )
    );
}

Foam::tmp<Foam::volSymmTensorField>
Foam::RASModels::SATFMcontinuousModel::devRhoReff() const
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
          + 2.0 * alpha_ * rho_ *
            symm(
                     (k_&eX)*(eX*eX)
                   + (k_&eY)*(eY*eY)
                   + (k_&eZ)*(eZ*eZ)
                 )
        )
    );
}


Foam::tmp<Foam::fvVectorMatrix>
Foam::RASModels::SATFMcontinuousModel::divDevRhoReff
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
    
    return
    (
      - fvm::laplacian(rho_*nut_, U)
      - fvc::div
        (
            (rho_*nut_)*dev2(T(fvc::grad(U)))
          - 2.0 * alpha_ * rho_ *
            (
                (k_&eX)*(eX*eX)
              + (k_&eY)*(eY*eY)
              + (k_&eZ)*(eZ*eZ)
            )
        )
    );
}

void Foam::RASModels::SATFMcontinuousModel::boundNormalStress
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

void Foam::RASModels::SATFMcontinuousModel::boundxiPhiG
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

void Foam::RASModels::SATFMcontinuousModel::boundS
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



void Foam::RASModels::SATFMcontinuousModel::correct()
{
    // Local references
    const twoPhaseSystem& fluid = refCast<const twoPhaseSystem>(phase_.fluid());
    volScalarField alpha(max(alpha_, scalar(0)));
    // solid volume fraction
    volScalarField alpha1 = 1.0 - alpha;
    const volScalarField& rho = phase_.rho();
    const volScalarField& rho1 = fluid.otherPhase(phase_).rho();
    const surfaceScalarField& alphaRhoPhi = alphaRhoPhi_;
    const volVectorField& U = U_;
    
    // dispersed Phase velocity
    const volVectorField& Ud_ = fluid.otherPhase(phase_).U();
    
    // slip velocity
    volVectorField uSlip = U - Ud_;
    
    // gravity vector
    const uniformDimensionedVectorField& g = mesh_.lookupObject<uniformDimensionedVectorField>("g");
    
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
    // gradient of continuous phase volume fraction
    volVectorField gradAlpha  = fvc::grad(alpha);
    
    // get turbulent kinetic energy of continuous-phase
    const volVectorField& kD_(mesh_.lookupObject<volVectorField>
                                     ("k." + fluid.otherPhase(phase_).name()));
    
    // get alphaP2Mean
    const volScalarField& alphaP2Mean1_(mesh_.lookupObject<volScalarField>
                             ("alphaP2Mean." + fluid.otherPhase(phase_).name()));
    volScalarField alphaP2MeanO = max(alphaP2Mean1_,alphaP2Mean_);
    
    const cellList& cells = mesh_.cells();
    
    // simple filter for smoothing of correlation coefficients
    simpleFilter filterS(mesh_);
    
    // get drag coefficient
    volScalarField beta
    (
        fluid.lookupSubModel<dragModel>
        (
            fluid.otherPhase(phase_),
            phase_
        ).K()
    );
    beta.max(1.0e-7);
    
    volScalarField betaA = beta/(rho*alpha);
    
    // get drift velocity
    volVectorField KdUdrift
    (
        fluid.lookupSubModel<driftVelocityModel>
        (
            fluid.otherPhase(phase_),
            phase_
        ).KdUdrift()
    );
    
    // compute total k
    volScalarField km  = k_ & eSum;
    km.max(kSmall.value());
    
    // compute grid size for mixing length
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
        deltaF_[cellI] = 2.0*deltaMaxTmp;
    }
    
    volScalarField wD = wallDist(mesh_).y();
    
    // correction for cases w/o walls
    // (since wall distance is then negative)
    deltaF_ = neg(wD)*deltaF_ + pos(wD)*min(deltaF_,2.0*wD);
    
    if (dynamicAdjustment_) {
        // precompute \bar phi
        volScalarField alpha2f = filter_(alpha);
        alpha2f.min(1.0);
        volScalarField alpha1f = scalar(1.0) - alpha2f;
        alpha1f.max(1.e-7);
        volScalarField alpha1fP2 = filter_(sqr(alpha1));
        alpha1fP2.max(sqr(residualAlpha_.value()));
        
        volVectorField Uf = filter_(alpha*U)/alpha2f;
        // compute xiPhiG_
        xiPhiG_ = - filterS(
                      filter_(alpha*U)
                    - alpha2f*filter_(U)
                  );
        volVectorField xiPhiGDenomSqr =
                  filterS(
                      (alpha1fP2-sqr(alpha1f))
                    * (
                          filter_(
                                     magSqr(U&eX)*eX
                                   + magSqr(U&eY)*eY
                                   + magSqr(U&eZ)*eZ
                                 )
                               / alpha2f
                        - magSqr(Uf&eX)*eX
                        - magSqr(Uf&eY)*eY
                        - magSqr(Uf&eZ)*eZ
                       )
                      );
        forAll(cells,cellI)
        {
            for (int i=0; i< 3; i++) {
                xiPhiG_[cellI].component(i) /= Foam::sqrt(Foam::max(xiPhiGDenomSqr[cellI].component(i),1.0e-14));
            }
        }
        // compute triple correlation
        volVectorField Ucf = filter_(alpha*U)/alpha2f;
        xiPhiGG_ = filterS(
                       filter_(alpha1*(U&U))
                     - alpha1f*filter_(U&U)
                     - 2.0*(Ucf&(filter_(alpha1*U) - alpha1f*filter_(U)))
                   )
                 / filterS(
                       sqrt(max(alpha1fP2-sqr(alpha1f),sqr(residualAlpha_)))
                     * max(mag(filter_(alpha*U&U)/alpha2f)-magSqr(Ucf),sqr(uSmall))
                   );

        // compute correlation coefficient between gas phase and solid phase velocity
        xiGS_ = filterS(
                     filter_(alpha1*(U&Ud_))/alpha1f
                   - (filter_(alpha1*U) & filter_(alpha1*Ud_))/sqr(alpha1f)
                )
              / filterS(
                    sqrt(max(filter_(alpha1*(U&U))/alpha1f-magSqr(filter_(alpha1*U)/alpha1f),kSmall))
                  * sqrt(max(filter_(alpha1*(Ud_&Ud_))/alpha1f-magSqr(filter_(alpha1*Ud_)/alpha1f),kSmall))
                 );

        // limit and smooth correlation coefficients
        // xiPhiG_
        xiPhiG_ = 0.5*(-mag(xiPhiG_)*gradAlpha
                       /(mag(gradAlpha)+dimensionedScalar("small",dimensionSet(0,-1,0,0,0),1.e-7)) + xiPhiG_);
        // xiPhiG_ = 0.5*(-mag(xiPhiG_)*uSlip/(mag(uSlip)+uSmall) + xiPhiG_);
        boundxiPhiG(xiPhiG_);
        
        // xiPhiGG_
        xiPhiGG_ = 0.5*(mag(xiPhiGG_) + xiPhiGG_);
        xiPhiGG_.max(-0.99);
        xiPhiGG_.min(0.99);
        
        // xiGS_ (xiGS_ is positive)
        xiGS_ = 0.5*(mag(xiGS_) + xiGS_);
        xiGS_.max(-1.0);
        xiGS_.min(1.0);
        
        // Currently no dynamic procedure for Cmu and Ceps
        // Set Cmu
        Cmu_ = CmuScalar_;
        // Set Ceps
        Ceps_ = CepsScalar_;
        // Set Cp
        Cp_     = CpScalar_;
        
        // compute mixing length dynamically
        lm_ = sqrt(km/(sqrt(SijSij&&SijSij)+dimensionedScalar("small",dimensionSet(0,0,-2,0,0),1.e-7)));
        lm_ = min(Cmu_*deltaF_,lm_);
    } else {
        // the sign of xiPhiG should be opposite to the slip velocity
        xiPhiG_ =   xiPhiContScalar_
                  * (sign(uSlip&eX)*eX + sign(uSlip&eY)*eY + sign(uSlip&eZ)*eZ)
                  * alpha;
        xiPhiGG_ = scalar(0.0);
        xiGS_   = xiGSScalar_;
        Cmu_    = CmuScalar_;
        Ceps_   = CepsScalar_;
        Cp_     = CpScalar_;
        
        // compute mixing length
        lm_ = Cmu_*deltaF_;
    }
    // contrain mixing length
    lm_.max(lSmall.value());
    // compute xiGatS
    xiGatS_ =  scalar(1.0) + xiPhiGG_*sqrt(alphaP2MeanO)
            / max(alpha1*alpha*(scalar(1.0) - xiPhiGG_*sqrt(alphaP2MeanO)/alpha),residualAlpha_);
    // xiGatS_ = filterS(xiGatS_);
    xiGatS_.max(1.0e-7);
    xiGatS_.min(2.0);
    // correct xiGS_
    xiGS_ *= sqrt(xiGatS_);

    // Compute k_
    // ---------------------------
    volVectorField pDil = Cp_*sqr(alpha)*alpha1*(rho1-rho)*g/beta;
    if (!equilibrium_) {
        fv::options& fvOptions(fv::options::New(mesh_));

        // Construct the transport equation for k
        // --> Stefanie
        Info << "Solving k-equation (continuous phase) ... " << endl;
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
                                 / (2.0 * sigma_)
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
                        sqrt((kD_&eX)*(k_&eX))*eX
                      + sqrt((kD_&eY)*(k_&eY))*eY
                      + sqrt((kD_&eZ)*(k_&eZ))*eZ
                    )
                )
          + fvm::Sp(-2.0*beta*xiGatS_,k_)
          // drag production and pressure dilation
          - (KdUdrift&eX)*((uSlip&eX) + (pDil&eX))*eX
          - (KdUdrift&eY)*((uSlip&eY) + (pDil&eY))*eY
          - (KdUdrift&eZ)*((uSlip&eZ) + (pDil&eZ))*eZ
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
        Ceps_ = CepsScalar_;
        // Equilibrium => dissipation == production
        // Schneiderbauer (2017), equ. (56)
        forAll(cells,cellI)
        {
            for (int i=0; i<3; i++) {
                k_[cellI].component(i) =
                    sqr(
                         - xiGatS_[cellI]*betaA[cellI]*lm_[cellI]
                         + Foam::sqrt(
                              sqr(xiGatS_[cellI]*betaA[cellI]*lm_[cellI])
                            + 2.0 * lm_[cellI]
                            * Foam::max(
                                 lm_[cellI]*SijSijV[cellI].component(i)
                               + betaA[cellI]*xiGS_[cellI]*Foam::sqrt(Foam::max(kD_[cellI].component(i),kSmall.value()))
                               - KdUdrift[cellI].component(i)*(uSlip[cellI].component(i) + pDil[cellI].component(i))
                                        /(2.0*alpha[cellI]*rho[cellI]*Foam::sqrt(Foam::max(k_[cellI].component(i),kSmall.value())))
                                , 0.
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
    //--------------------------------------------
    // update km
    km = k_ & eSum;
    km.max(kSmall.value());
    volScalarField divU(fvc::div(U));
    volScalarField denom = divU + CphiGscalar_ * Ceps_ * sqrt(km)/lm_;
    denom.max(kSmall.value());
    
    Info << "Computing alphaP2Mean (continuous phase) ... " << endl;
    if (dynamicAdjustment_) {
        volVectorField xiKgradAlpha = (
                                         ((sqrt(k_&eX) * (gradAlpha&eX) * (xiPhiG_&eX)) * eX)
                                       + ((sqrt(k_&eY) * (gradAlpha&eY) * (xiPhiG_&eY)) * eY)
                                       + ((sqrt(k_&eZ) * (gradAlpha&eZ) * (xiPhiG_&eZ)) * eZ)
                                       );
        alphaP2Mean_ =   8.0
                       * magSqr(xiKgradAlpha)
                       / sqr(denom);
    } else {
        alphaP2Mean_ =   8.0
                       * magSqr(xiPhiG_)
                       * sqr(
                                (sqrt(k_&eX) * mag(gradAlpha&eX))
                              + (sqrt(k_&eY) * mag(gradAlpha&eY))
                              + (sqrt(k_&eZ) * mag(gradAlpha&eZ))
                         )
                       / sqr(denom);
    }
    // limti alphaP2Mean_
    alphaP2Mean_.max(sqr(residualAlpha_.value()));
    alphaP2Mean_ = min(alphaP2Mean_, alpha*(1.0 - alpha));
    
    // compute nut_ (Schneiderbauer, 2017; equ. (34))
    nut_ = pos((scalar(1.0) - alpha) - residualAlpha_)*alpha*sqrt(km)*lm_;
    
    // Limit viscosity and add frictional viscosity
    nut_.min(maxNut_);
    
    Info<< "SA-TFM (continuous Phase):" << nl
    << "    max(nut) = " << max(nut_).value() << nl
    << "    max(k_)  = " << max(km).value() << endl;
}


// ************************************************************************* //

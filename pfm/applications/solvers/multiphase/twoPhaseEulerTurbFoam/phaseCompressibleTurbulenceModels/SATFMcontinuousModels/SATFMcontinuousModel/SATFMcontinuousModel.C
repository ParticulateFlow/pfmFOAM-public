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
    anIsoTropicNut_(coeffDict_.lookup("anIsoTropicNut")),

    alphaMax_("alphaMax", dimless, coeffDict_),

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
        coeffDict_.lookupOrDefault<scalar>("xiPhiG",0.5)
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
        coeffDict_.lookupOrDefault<scalar>("CphiG",0.2)
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

    xiPhiG_
    (
        IOobject
        (
            "xiPhiG",
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedVector("value", dimensionSet(0, 0, 0, 0, 0), vector(-0.5,-0.5,-0.5)),
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
            IOobject::AUTO_WRITE
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
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("value", dimensionSet(0, 0, 0, 0, 0), 0.4),
        // Set Boundary condition
        zeroGradientFvPatchField<scalar>::typeName
    ),

    CphiG_
    (
        IOobject
        (
            IOobject::groupName("CphiG", phase.name()),
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("value", dimensionSet(0, 0, 0, 0, 0), 0.2),
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


    R2_
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
        coeffDict().lookup("anIsoTropicNut") >> anIsoTropicNut_;
        alphaMax_.readIfPresent(coeffDict());
        xiPhiContScalar_.readIfPresent(coeffDict());
        CmuScalar_.readIfPresent(coeffDict());
        CphiGscalar_.readIfPresent(coeffDict());
        CepsScalar_.readIfPresent(coeffDict());
        CpScalar_.readIfPresent(coeffDict());
        sigma_.readIfPresent(coeffDict());
        gN_.readIfPresent(coeffDict());
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
          + 2.0 * alpha_ * symm(R2_)
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
          + 2.0 * alpha_ * symm(R2_)
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
    pos(scalar(1.0) - alpha_ - residualAlpha_)*
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
                 R2_
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

void Foam::RASModels::SATFMcontinuousModel::boundCorrTensor
(
    volTensorField& R
) const
{
    scalar xiMin = -1.0;
    scalar xiMax = 1.0;

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

void Foam::RASModels::SATFMcontinuousModel::boundGradU
(
    volTensorField& R
) const
{
    scalar sMin = -1.0e2;
    scalar sMax =  1.0e2;

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
    volScalarField alpha(max(alpha_, residualAlpha_));
    // solid volume fraction
    volScalarField alpha1 = max(1.0 - alpha,residualAlpha_);
    const volScalarField& rho = phase_.rho();
    const volScalarField& rho1 = fluid.otherPhase(phase_).rho();
    const surfaceScalarField& alphaRhoPhi = alphaRhoPhi_;
    const surfaceScalarField& phi2 = phi_;
    const volVectorField& U = U_;
    
    // dispersed Phase velocity
    const volVectorField& Ud_ = fluid.otherPhase(phase_).U();
    
    // slip velocity
    volVectorField uSlip = U - Ud_;
    
    // gravity vector
    //const uniformDimensionedVectorField& g = mesh_.lookupObject<uniformDimensionedVectorField>("g");
    
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

    // gradient of continuous phase volume fraction
    volVectorField gradAlpha  = fvc::grad(alpha);
    
    // get turbulent kinetic energy of continuous-phase
    const volVectorField& kD_(mesh_.lookupObject<volVectorField>
                                     ("k." + fluid.otherPhase(phase_).name()));
    
    // get correlation coefficient between continuous and solid phase velocities
    const volScalarField& xiGS_(mesh_.lookupObject<volScalarField>
                             ("xiGS"));
    // get correction coefficient for k seen by the particles
    const volScalarField& xiGatS_(mesh_.lookupObject<volScalarField>
                             ("xiGatS"));
    
    const cellList& cells = mesh_.cells();
    
    // simple filter for local smoothing
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
    volScalarField km  = (k_ & eSum);
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
    deltaF_.max(lSmall.value());
    
    if (dynamicAdjustment_) {
        // precompute \bar phi
        volScalarField alpha2f = filter_(alpha);
        alpha2f.min(1.0);
        volScalarField alpha1f = scalar(1.0) - alpha2f;
        alpha1f.max(residualAlpha_.value());
        volScalarField alpha1fP2 = filter_(sqr(alpha1));
        alpha1fP2.max(sqr(residualAlpha_.value()));
        
        volVectorField Uf = filter_(alpha*U)/alpha2f;
        // compute xiPhiG_
        volVectorField xiPhiGNom =   (
                                          filter_(alpha1*U)
                                        - alpha1f*filter_(U)
                                      );
        volScalarField tmpA = alpha1fP2-sqr(alpha1f);
        volScalarField tmpDenX = tmpA
                              * (
                                    filter_(alpha*sqr(U&eX)) / alpha2f
                                  - sqr(Uf&eX)
                                );
        volScalarField tmpDenY = tmpA
                               * (
                                    filter_(alpha*sqr(U&eY)) / alpha2f
                                  - sqr(Uf&eY)
                                 );
        volScalarField tmpDenZ = tmpA
                               * (
                                    filter_(alpha*sqr(U&eZ)) / alpha2f
                                  - sqr(Uf&eZ)
                                );
       
        tmpDenX.max(ROOTVSMALL);
        tmpDenY.max(ROOTVSMALL);
        tmpDenZ.max(ROOTVSMALL);

        xiPhiG_ =  eX
                 * (
                        filterS(sqrt(tmpDenX)*(xiPhiGNom&eX))/filterS(tmpDenX)
                    )
                 + eY
                 * (
                        filterS(sqrt(tmpDenY)*(xiPhiGNom&eY))/filterS(tmpDenY)
                    )
                 + eZ
                 * (
                        filterS(sqrt(tmpDenZ)*(xiPhiGNom&eZ))/filterS(tmpDenZ)
                    );
        /*
        volScalarField xiPhiGDenomSqr =   (alpha1fP2-sqr(alpha1f))
                                        * (filter_(alpha*magSqr(U))/alpha2f - magSqr(Uf));
        xiPhiGDenomSqr.max(ROOTVSMALL);
        xiPhiG_ = 3.0*filterS(xiPhiGNom*sqrt(xiPhiGDenomSqr))/filterS(xiPhiGDenomSqr);
        */
        // limit and smooth correlation coefficients
        // xiPhiG_
        /*
        xiPhiG_ = 0.5*(
                        - mag(xiPhiG_)*gradAlpha
                              /(mag(gradAlpha)+dimensionedScalar("small",dimensionSet(0,-1,0,0,0),1.e-7))
                        + xiPhiG_
                      );
        */
        boundxiPhiG(xiPhiG_);
               
        // compute mixing length dynamically
        /*
        volScalarField Lij      = filter_(alpha*magSqr(U))/alpha2f - magSqr(Uf);
        volScalarField magSqrDf = filter_(alpha*magSqr(D))/alpha2f;
        magSqrDf.max(SMALL);
        volSymmTensorField Df   = filter_(alpha*D)/alpha2f;
        volScalarField Mij      = sqr(deltaF_)*(4.0*magSqr(Df) - magSqrDf);
        volScalarField MijMij   = filterS(sqr(Mij));
        MijMij.max(SMALL);
        
        volScalarField CmuT     = 0.5*(filterS(Lij * Mij)/(MijMij));
        CmuT.min(2.0*sqr(CmuScalar_.value()));
        CmuT.max(0.01);
        
        Cmu_ = pos(scalar(1.0) - alpha_ - residualAlpha_)*sqrt(CmuT)
             + neg(scalar(1.0) - alpha_ - residualAlpha_)*CmuScalar_;
        */
        Cmu_    = CmuScalar_;
        // dynamic procedure for Ceps
        /*
        volScalarField nu2 = mesh_.lookupObject<volScalarField>("thermo:mu." + phase_.name())/rho_;
        volScalarField magSqrD = magSqr(D);
        volScalarField LijEps = nu2*alpha2f*(magSqrDf - magSqr(Df));
        volScalarField MijEps = sqr(deltaF_*Cmu_)
                              * (
                                    alpha2f*magSqrDf*sqrt(magSqrDf)
                                  - filter_(alpha*magSqrD*sqrt(magSqrD))
                                );
        volScalarField MijMijEps = filterS(sqr(MijEps));
        MijMijEps.max(SMALL);
        
        volScalarField CepsT = 2.0*mag(filterS(LijEps * MijEps)/(MijMijEps));
        
        Ceps_ = pos(scalar(1.0) - alpha_ - residualAlpha_)*(CepsT)
              + neg(scalar(1.0) - alpha_ - residualAlpha_);
        
        Ceps_.min(1.0);
        Ceps_.max(0.01);
        */
        Ceps_   = CepsScalar_;
        // Compute CphiG_
        CphiG_ = CphiGscalar_/Ceps_;
        
        // Currently no dynamic procedure for Cp
        /*
        const volScalarField& p(mesh_.lookupObject<volScalarField>("p"));
        volScalarField rhom = rho*alpha + alpha1*rho1;
        volVectorField gradp = fvc::grad(p);

        Cp_ = (gradp&gN_)/(rhom*(gN_&gN_));
        Cp_ = filterS(Cp_);
        Cp_.min(1.0);
        Cp_.max(0.1);
        */
        Cp_     = CpScalar_;
    } else {
        // the sign of xiPhiG should be opposite to the slip velocity
        volVectorField xiPhiGDir = uSlip/(mag(uSlip)+dimensionedScalar("small",dimensionSet(0,1,-1,0,0),1.e-7));
        xiPhiG_ = - (xiPhiContScalar_)
                  * (scalar(1.0) - alpha1)
                  * xiPhiGDir;
        Cmu_    = CmuScalar_;
        Ceps_   = CepsScalar_;
        Cp_     = CpScalar_;
        CphiG_  = CphiGscalar_/Ceps_;
    }
    // compute mixing length
    lm_ = Cmu_*deltaF_;
    


    // Compute k_
    // ---------------------------
    volVectorField pDil = Cp_*sqr(alpha)*alpha1*(rho1-rho)*gN_/beta;
    if (!equilibrium_) {
        if (anIsoTropicNut_) {
            volTensorField gradUR2 = 0.5*((R2_&gradU) + ((R2_.T())&(gradU.T())));
            volVectorField shearProd_ =   (gradUR2&&(eX*eX))*(eX)
                                        + (gradUR2&&(eY*eY))*(eY)
                                        + (gradUR2&&(eZ*eZ))*(eZ);
        } else {
            volVectorField shearProd_ = -lm_*(
                                                 (((SijSij&eX)&eSum)*sqrt(k_&eX))*eX
                                               + (((SijSij&eY)&eSum)*sqrt(k_&eY))*eY
                                               + (((SijSij&eZ)&eSum)*sqrt(k_&eZ))*eZ
                                             );
        }
        
        fv::options& fvOptions(fv::options::New(mesh_));

        // Construct the transport equation for k
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
                                 / (sigma_)
                           , k_, "laplacian(kappa,k)")
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
                               - KdUdrift[cellI].component(i)*(uSlip[cellI].component(i))
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
    km = (k_ & eSum);
    km.max(kSmall.value());

    // compute nut_
    nut_ = alpha*sqrt(km)*lm_;
    
    // compute fields for transport equation for phiP2
    volScalarField divU(fvc::div(U));
    volScalarField dissPhiP2 = CphiG_ * Ceps_ * sqrt(km)/lm_;
    volScalarField denom = mag(divU) + dissPhiP2;
    denom.max(SMALL);
    volScalarField xiKgradAlpha = - (
                                       ((sqrt(k_&eX) * (gradAlpha&eX) * (xiPhiG_&eX)))
                                     + ((sqrt(k_&eY) * (gradAlpha&eY) * (xiPhiG_&eY)))
                                     + ((sqrt(k_&eZ) * (gradAlpha&eZ) * (xiPhiG_&eZ)))
                                    );
    
    Info << "Computing alphaP2Mean (continuous phase) ... " << endl;
    if (!equilibrium_) {
        // Construct the transport equation for alphaP2Mean
        fvScalarMatrix phiP2Eqn
        (
            fvm::ddt(alphaP2Mean_)
          + fvm::div(phi2, alphaP2Mean_)
          // diffusion
          - fvc::div(
                        (
                           alpha1*lm_
                         * (
                              (sqrt(k_&eX)*(eX*eX))
                            + (sqrt(k_&eY)*(eY*eY))
                            + (sqrt(k_&eZ)*(eZ*eZ))
                           )
                         / (sigma_)
                       )
                     & (fvc::grad(alphaP2Mean_/alpha1))
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
          //+ 2.0*nut_*magSqr(gradAlpha)
          // production/dissipation
          + fvm::Sp(-dissPhiP2,alphaP2Mean_)
        );

        phiP2Eqn.relax();
        phiP2Eqn.solve();
    } else {
        alphaP2Mean_ =   8.0
                       * sqr(xiKgradAlpha)
                       / sqr(denom)
                       * neg(xiKgradAlpha);
    }
    // limti alphaP2Mean_
    volScalarField alphaM  = alphaMax_ - alpha1;
    alphaM.max(0.0);
    volScalarField alphaL2 = sqr(min(alpha1,alphaM));
    alphaP2Mean_ = min(
                         alphaP2Mean_,
                         0.99*alphaL2
                      );
    alphaP2Mean_.max(ROOTVSMALL);
    alphaP2Mean_.correctBoundaryConditions();

    
    if (anIsoTropicNut_) {
        nut_ = dimensionedScalar("zero", dimensionSet(0, 2, -1, 0, 0), 0.0);
        volScalarField alphaf = filter_(alpha);
        alphaf.max(residualAlpha_.value());
        volVectorField Uf = filter_(alpha*U)/alphaf;
        // compute correlation coefficients
        xiUU_ = (filter_(alpha*(U*U))/alphaf - Uf*Uf)
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
                        R2_[cellI].component(j+i*3) =  (xiUU_[cellI].component(j+i*3))
                                *sqrt(k_[cellI].component(i)*k_[cellI].component(j));
                    } else {
                        R2_[cellI].component(j+i*3) =  sqrt(k_[cellI].component(i)*k_[cellI].component(j));
                    }
                }
            }
        }
    } else {
        // compute R2_
        R2_ = zeroR;
    }
    
    // Limit viscosity and add frictional viscosity
    nut_.min(maxNut_);
    
    Info << "SA-TFM (continuous Phase):" << nl
         << "    max(nut)        = " << max(nut_).value() << nl
         << "    max(phiP2/phi2) = " << max(alphaP2Mean_/sqr(alpha1)).value() << nl
         << "    max(k2)         = " << max(k_&eSum).value() << nl
         << "    mean(k2x)       = " << fvc::domainIntegrate(alpha*(k_&eX)).value()
                                        /fvc::domainIntegrate(alpha).value()
                                     << nl
         << "    mean(k2y)       = " << fvc::domainIntegrate(alpha*(k_&eY)).value()
                                        /fvc::domainIntegrate(alpha).value()
                                     << nl
         << "    mean(k2z)       = " << fvc::domainIntegrate(alpha*(k_&eZ)).value()
                                        /fvc::domainIntegrate(alpha).value()
                                     << nl    
         << "    mean(phi2P2)    = " << fvc::domainIntegrate(alpha*(alphaP2Mean_)).value()
                                        /fvc::domainIntegrate(alpha).value()
         << endl;
}


// ************************************************************************* //

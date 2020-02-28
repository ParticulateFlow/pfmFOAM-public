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

    equilibriumK_(coeffDict_.lookup("equilibriumK")),
    equilibriumPhiP2_(coeffDict_.lookup("equilibriumPhiP2")),
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
        coeffDict_.lookupOrDefault<scalar>("maxNut",1)
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
        dimensionedScalar("value", dimensionSet(0, 0, 0, 0, 0), 0.4)
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
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("value", dimensionSet(0, 0, 0, 0, 0), 0.4)
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
        dimensionedScalar("value", dimensionSet(0, 0, 0, 0, 0), 0.2)
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

    R2_
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
        coeffDict().lookup("equilibriumK") >> equilibriumK_;
        coeffDict().lookup("equilibriumPhiP2") >> equilibriumPhiP2_;
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
    dimensionedScalar uSmall("uSmall", U_.dimensions(), 1.0e-6);
    dimensionedScalar kSmall("kSmall", k_.dimensions(), 1.0e-6);
    
    tmp<volScalarField> kT
    (
        
       1.5
      *min(
           max
          (
                (k_&eSum)
              - mag(k_ & U_)
              / (mag(U_)+uSmall)
            ,
                kSmall
           )
      ,3.0*maxK_)
     );
    return kT;
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::SATFMcontinuousModel::epsilon() const
{
    dimensionedVector eSum
    (
        "eSum",
        dimensionSet(0, 0, 0, 0, 0, 0, 0),
        vector(1,1,1)
    );
    return Ceps_*pow(k(),3.0/2.0)/lm_;
}


Foam::tmp<Foam::volSymmTensorField>
Foam::RASModels::SATFMcontinuousModel::R() const
{
    if (anIsoTropicNut_) {
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
                2.0 * alpha_ * symm(R2_)
              - (nuEff()-nut_)*(twoSymm(fvc::grad(U_)))
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
                2.0 * alpha_ * symm(R2_)
              - (nuEff())*(twoSymm(fvc::grad(U_)))
            )
        );
    }
}

Foam::tmp<Foam::volSymmTensorField>
Foam::RASModels::SATFMcontinuousModel::devRhoReff() const
{
    if (anIsoTropicNut_) {
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
                2.0 * alpha_ * rho_* symm(R2_)
              - rho_*(nuEff()-nut_)*dev(twoSymm(fvc::grad(U_)))
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
                2.0 * alpha_ * rho_ * symm(R2_)
              - rho_*(nuEff())*dev(twoSymm(fvc::grad(U_)))
            )
        );
    }
}


Foam::tmp<Foam::fvVectorMatrix>
Foam::RASModels::SATFMcontinuousModel::divDevRhoReff
(
    volVectorField& U
) const
{
    if (!anIsoTropicNut_) {
        return
        (
          - fvm::laplacian(rho_*nuEff(), U)
          - fvc::div
            (
                (rho_*nuEff())*dev2(T(fvc::grad(U)))
            )
         + fvc::div
           (
                2.0
              * alpha_
              * rho_
              * R2_
           )
        );
    } else {
        return
        (
          - fvm::laplacian(rho_*(nuEff()), U)
          - fvc::div
            (
                rho_*(nuEff() - nut_)*dev2(T(fvc::grad(U)))
            )
          + fvc::div
            (
                 2.0
               * alpha_
               * rho_
               * R2_
            )
          + fvc::laplacian(rho_*nut_, U)
        );
    }
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

void Foam::RASModels::SATFMcontinuousModel::boundCorrTensor
(
    volTensorField& R
) const
{
    scalar xiMin = -0.75;
    scalar xiMax = 0.75;

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




void Foam::RASModels::SATFMcontinuousModel::correct()
{
    // Local references
    const twoPhaseSystem& fluid = refCast<const twoPhaseSystem>(phase_.fluid());
    volScalarField alpha(min(max(alpha_, 1.e-7),1.0));
    // solid volume fraction
    volScalarField alpha1 = max(1.0 - alpha,1.e-7);
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
    volScalarField km(k_&eSum);
    km.max(kSmall.value());
    
    // local reference to deltaF
    volScalarField deltaF(deltaF_);
    // compute grid size for mixing length
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
        deltaF[cellI] = 2.0*deltaMaxTmp;
    }
    
    volScalarField wD = wallDist(mesh_).y();
    // correction for cases w/o walls
    // (since wall distance is then negative)
    deltaF_ = neg(wD)*deltaF + pos(wD)*min(deltaF,wD);
    // correction for cyclic patches
    {
        const fvPatchList& patches = mesh_.boundary();
        
        forAll(patches, patchi) {
            const fvPatch& curPatch = patches[patchi];

            if (isA<cyclicAMIFvPatch>(curPatch)||isA<cyclicACMIFvPatch>(curPatch)) {
                forAll(curPatch, facei) {
                    label celli = curPatch.faceCells()[facei];
                    deltaF_[celli] = deltaF[celli];
                }
            }
        }
    }
    deltaF_.max(lSmall.value());
    
    // compute nut
    nut_ = alpha*sqrt(km)*lm_;
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
        tmpA.max(ROOTVSMALL);

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
        /*
        xiPhiG_ =  eX
                 * (
                        ((xiPhiGNom&eX))/sqrt(tmpDenX)
                    )
                 + eY
                 * (
                        ((xiPhiGNom&eY))/sqrt(tmpDenY)
                    )
                 + eZ
                 * (
                        ((xiPhiGNom&eZ))/sqrt(tmpDenZ)
                    );
         */
        xiPhiG_ =  eX
                 * (
                        filterS((xiPhiGNom&eX)*sqrt(tmpDenX))/filterS(tmpDenX)
                    )
                 + eY
                 * (
                        filterS((xiPhiGNom&eY)*sqrt(tmpDenY))/filterS(tmpDenY)
                    )
                 + eZ
                 * (
                        filterS((xiPhiGNom&eZ)*sqrt(tmpDenZ))/filterS(tmpDenZ)
                    );
        // wall treatment for xiPhiG
        /*
        const fvPatchList& patches = mesh_.boundary();
        volVectorField::Boundary& xiPhiGBf = xiPhiG_.boundaryFieldRef();
        
        forAll(patches, patchi) {
            const fvPatch& curPatch = patches[patchi];

            if (isA<wallFvPatch>(curPatch)) {
                vectorField& xiPhiGw = xiPhiGBf[patchi];
                forAll(curPatch, facei) {
                    label celli = curPatch.faceCells()[facei];
                    xiPhiGw[facei] = -xiPhiContScalar_.value()
                                     *uSlip[celli]
                                     /(mag(uSlip[celli])+SMALL);
                }
            }
        }
        */
        // xiPhiG_ = filterS(xiPhiG_);
        // limit xiPhiG_
        boundxiPhiG(xiPhiG_);
        
        // compute mixing length dynamically
        /*
        volScalarField Lij      = filter_(alpha*magSqr(U))/alpha2f - magSqr(Uf);
        Lij.max(SMALL);
        volScalarField magSqrDf = filter_(magSqr(D));
        magSqrDf.max(SMALL);
        volSymmTensorField Df   = filter_(D);
        volScalarField Mij      = sqr(deltaF_)*(4.0*magSqr(Df) - magSqrDf);
        volScalarField MijMij   = filterS(sqr(Mij));
        MijMij.max(SMALL);
        
        volScalarField CmuT     = 0.5*(filterS(Lij * Mij)/(MijMij));
        CmuT.min(4.0*sqr(CmuScalar_.value()));
        Cmu_ = sqrt(0.5*(CmuT+mag(CmuT))+scalar(1.e-2));
        */
        Cmu_    = CmuScalar_;
        // dynamic procedure for Ceps
        /*
        volScalarField LijEps    = alpha*nuEff()*(magSqrDf - magSqr(Df));
        volScalarField MijEps    = pow(alpha*Lij,1.5)/(2.0*lm_);
        volScalarField MijMijEps = filterS(sqr(MijEps));
        MijMijEps.max(SMALL);
        
        volScalarField CepsT     = filterS(LijEps*MijEps)/(MijMijEps);
        
        Ceps_ = 0.33*(CepsT + mag(CepsT) + CepsScalar_);
        
        Ceps_.min(2.0);
        */
        Ceps_ = CepsScalar_;
        // Compute CphiG_
        CphiG_ = CphiGscalar_;
        
        // Currently no dynamic procedure for Cp
        /*
        const volScalarField& p_rgh(mesh_.lookupObject<volScalarField>("p_rgh"));
        volScalarField rhom = rho*alpha + alpha1*rho1;
        volVectorField gradp = fvc::grad(p_rgh);

        Cp_ = filterS((gradp&gN_)*rhom)/filterS(sqr(rhom)*magSqr(gN_))
            + dimensionedScalar("unity",dimensionSet(0,0,0,0,0),1.0);
        Cp_ = 0.33*(Cp_ + mag(Cp_) + CpScalar_);
        Cp_.min(1.0);
        Cp_ = pos(scalar(1.0) - alpha_ - residualAlpha_)*Cp_
            + neg(scalar(1.0) - alpha_ - residualAlpha_);
        */
        Cp_ = CpScalar_;
    } else {
        // the sign of xiPhiG should be opposite to the slip velocity
        volVectorField xiPhiGDir = uSlip/(mag(uSlip)+dimensionedScalar("small",dimensionSet(0,1,-1,0,0),1.e-7));
        xiPhiG_ = - (xiPhiContScalar_)
                  * (scalar(1.0) - alpha1)
                  * xiPhiGDir;
        Cmu_    = CmuScalar_;
        Ceps_   = CepsScalar_;
        Cp_     = CpScalar_;
        CphiG_  = CphiGscalar_;
    }
    
    // compute mixing length
    lm_ = Cmu_*deltaF_;
    


    // Compute k_
    // ---------------------------
    volVectorField pDil = Cp_*sqr(alpha)*alpha1*(rho1-rho)*gN_/beta;
    if (!equilibriumK_) {
        // compute production term according to Reynolds-stress model
        volTensorField R2t(R2_);
        if (!anIsoTropicNut_) {
            R2t -= 0.5*nut_*dev(gradU + gradU.T());
        }
        volTensorField gradUR2 = 0.5*((R2t&gradU) + ((gradU.T())&(R2t.T())));

        //volTensorField gradUR2 = 0.5*((R2_&gradU) + ((R2_.T())&(gradU.T())));
        /*
        shearProd_ = pos(mag(wD) - deltaF)
                    *(
                         (gradUR2&&(eX*eX))*(eX)
                       + (gradUR2&&(eY*eY))*(eY)
                       + (gradUR2&&(eZ*eZ))*(eZ)
                      )
                    // special treatment of P_k near walls
                   - neg(mag(wD) - deltaF)
                    *lm_
                    *(
                         (sqr(U&eX)*sqrt(k_&eX))*eX
                       + (sqr(U&eY)*sqrt(k_&eY))*eY
                       + (sqr(U&eZ)*sqrt(k_&eZ))*eZ
                      )
                     /sqr(deltaF_);
        */
        shearProd_ = (
                           (gradUR2&&(eX*eX))*(eX)
                         + (gradUR2&&(eY*eY))*(eY)
                         + (gradUR2&&(eZ*eZ))*(eZ)
                      );
        // compute prefactor for dissipation term
        // volScalarField coeffDissipation(Ceps_*alpha*rho/lm_);
        
        fv::options& fvOptions(fv::options::New(mesh_));

        // Construct the transport equation for k
        Info << "Solving k-equation (continuous phase) ... " << endl;
        fvVectorMatrix kEqn
        (
            fvm::ddt(alpha, rho, k_)
          + fvm::div(alphaRhoPhi, k_)
          //- fvc::Sp((fvc::ddt(alpha, rho) + fvc::div(alphaRhoPhi)), k_)
          // diffusion with anisotropic diffusivity
          /*
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
        */
          - fvm::laplacian(
                             alpha*rho*sqrt(k())*lm_/(sigma_),
                             k_,
                             "laplacian(kappa,k)"
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
                        sqrt((kD_&eX)*(k_&eX))*eX
                      + sqrt((kD_&eY)*(k_&eY))*eY
                      + sqrt((kD_&eZ)*(k_&eZ))*eZ
                    )
                )
          + fvm::Sp(
                     - 2.0
                      *beta
                      *xiGatS_
                    ,
                      k_
                   )
          // drag production and pressure dilation
          - (KdUdrift&eX)*((uSlip&eX) + (pDil&eX))*eX
          - (KdUdrift&eY)*((uSlip&eY) + (pDil&eY))*eY
          - (KdUdrift&eZ)*((uSlip&eZ) + (pDil&eZ))*eZ
          // dissipation
          // - coeffDissipation*(sqrt(k_&eX)*(k_&eX))*eX
          // - coeffDissipation*(sqrt(k_&eY)*(k_&eY))*eY
          // - coeffDissipation*(sqrt(k_&eZ)*(k_&eZ))*eZ
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
        volVectorField SijSijV =  pos(mag(wD) - deltaF)
                                 *(
                                      ((SijSij&eX)&eSum)*eX
                                    + ((SijSij&eY)&eSum)*eY
                                    + ((SijSij&eZ)&eSum)*eZ
                                   )
                                 // special treatment of P_k near walls
                                + neg(mag(wD) - deltaF)
                                 *(
                                      (sqr(U&eX))*eX
                                    + (sqr(U&eY))*eY
                                    + (sqr(U&eZ))*eZ
                                   )
                                  /sqr(deltaF_);

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
    // limit k
    boundNormalStress(k_);
    // correct BCs
    k_.correctBoundaryConditions();

    //- compute variance of solids volume fraction
    //--------------------------------------------
    // update km
    km = k_&eSum;
    km.max(kSmall.value());
    
    // compute fields for transport equation for phiP2
//    alphaP2Mean_ = mesh_.lookupObject<volScalarField>("alphaP2Mean." + fluid.otherPhase(phase_).name());
    volScalarField divU(fvc::div(U));
    volScalarField dissPhiP2 = CphiG_ * Ceps_ * sqrt(km)/lm_;
    volScalarField denom  = mag(divU) + dissPhiP2;
    denom.max(SMALL);
    volScalarField xiKgradAlpha = - (
                                       ((sqrt(k_&eX) * (gradAlpha&eX) * (xiPhiG_&eX)))
                                     + ((sqrt(k_&eY) * (gradAlpha&eY) * (xiPhiG_&eY)))
                                     + ((sqrt(k_&eZ) * (gradAlpha&eZ) * (xiPhiG_&eZ)))
                                    )
                                  + CphiG_
                                   *(
                                        (KdUdrift&uSlip)
                                      + (KdUdrift&pDil)
                                      - 2.0*beta*(xiGS_*sqrt(km*(kD_&eSum)) - xiGatS_*km)
                                    )
                                   /(rho*km);

    Info << "Computing alphaP2Mean (continuous phase) ... " << endl;
    if (!equilibriumPhiP2_) {
        // Construct the transport equation for alphaP2Mean
        fvScalarMatrix phiP2Eqn
        (
            fvm::ddt(alphaP2Mean_)
          + fvm::div(phi2, alphaP2Mean_)
          // diffusion
         /*
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
          */
          - fvc::div(
                       alpha*sqrt(k())*lm_/(sigma_)
                     * fvc::grad(alphaP2Mean_/alpha)
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
    alpha1.min(alphaMax_.value());
    volScalarField alphaM(alphaMax_-alpha1);
    alphaM.max(0);
    volScalarField alphaL2(alpha1*alphaM);
//    volScalarField cbrtPhiPhiM(cbrt(alpha1/alphaMax_));
//    volScalarField alphaL2 = sqr(alpha1)//*alphaM
//                            *(scalar(1.0) + cbrtPhiPhiM)
//                            /(
//                                 scalar(1.0)
//                               + (1.0/3.0)
//                                *cbrtPhiPhiM
//                                /(scalar(1.0) - cbrtPhiPhiM)
//                             );//sqr(min(alpha1,alphaM));
    alphaP2Mean_ = min(
                         alphaP2Mean_,
                         alphaL2
                      );
    alphaP2Mean_.max(VSMALL);
    alphaP2Mean_.correctBoundaryConditions();
    
    // use k() for nut in stress tensor
    nut_ = alpha*sqrt(k())*lm_;
    // Limit viscosity
    nut_.min(maxNut_);
    nut_.correctBoundaryConditions();
       
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
                        R2_[cellI].component(j+i*3) =  (xiUU_[cellI].component(j+i*3))
                                *sqrt(k_[cellI].component(i)*k_[cellI].component(j));
                    } else {
                        R2_[cellI].component(j+i*3) =  sqrt(k_[cellI].component(i)*k_[cellI].component(j));
                    }
                }
            }
        }
        
        // set wall-bc for R
        const fvPatchList& patches = mesh_.boundary();

        volTensorField::Boundary& RBf = R2_.boundaryFieldRef();

        forAll(patches, patchi) {
            const fvPatch& curPatch = patches[patchi];

            if (isA<wallFvPatch>(curPatch)) {
                tensorField& Rw = RBf[patchi];

                const scalarField& nutw = nut_.boundaryField()[patchi];

                const vectorField snGradU
                (
                    U.boundaryField()[patchi].snGrad()
                );

                const vectorField& faceAreas
                    = mesh_.Sf().boundaryField()[patchi];

                const scalarField& magFaceAreas
                    = mesh_.magSf().boundaryField()[patchi];

                forAll(curPatch, facei)
                {
                    // Calculate near-wall velocity gradient
                    const tensor gradUw
                        = (faceAreas[facei]/magFaceAreas[facei])*snGradU[facei];

                    // Set the wall Reynolds-stress to the near-wall shear-stress
                    // Note: the spherical part of the normal stress is included in
                    // the pressure
                    Rw[facei] = -nutw[facei]*(gradUw + gradUw.T());
                }
            }
        }
    } else {
        R2_  = 0*((k_&eX)*(eX*eX) + (k_&eY)*(eY*eY) + (k_&eZ)*(eZ*eZ));
    }
    
    R2_.correctBoundaryConditions();
    

    
    Info << "SA-TFM (continuous Phase):" << nl
         << "    max(nuEff)      = " << max(nuEff()).value() << nl
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
         << endl;
}


// ************************************************************************* //

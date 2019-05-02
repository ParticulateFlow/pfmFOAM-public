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

    CepsScalar_
    (
        "CepsScalar",
        dimensionSet(0,0,0,0,0),
        coeffDict_.lookupOrDefault<scalar>("Ceps",1.0)
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
        dimensionedVector("value", dimensionSet(0, 0, 0, 0, 0), vector(-0.5,-0.5,-0.5))
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
        dimensionedScalar("value", dimensionSet(0, 0, 0, 0, 0), 0.9)
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
        dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.0)
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
            IOobject::NO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("value", dimensionSet(0, 0, 0, 0, 0), 1.0)
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
        CepsScalar_.readIfPresent(coeffDict());
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
          + 2.0 * pos((scalar(1.0) - alpha_) - residualAlpha_) * alpha_ * rho_ *
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
          + 2.0 * pos((scalar(1.0) - alpha_) - residualAlpha_) * alpha_ * rho_ *
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
    pos((scalar(1.0) - alpha_) - residualAlpha_)*
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


void Foam::RASModels::SATFMcontinuousModel::correct()
{
    // Local references
    const twoPhaseSystem& fluid = refCast<const twoPhaseSystem>(phase_.fluid());
    volScalarField alpha(max(alpha_, scalar(0)));
    const volScalarField& rho = phase_.rho();
    const surfaceScalarField& alphaRhoPhi = alphaRhoPhi_;
    const volVectorField& U = U_;
    
    // dispersed Phase velocity
    const volVectorField& Ud_ = fluid.otherPhase(phase_).U();
    
    // slip velocity
    volVectorField uSlip = U - Ud_;
    
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
    volSymmTensorField D(symm(gradU));
    
    // compute S_{ij}S_{ij} (no summation over i!!)
    volVectorField SijSij =  magSqr(gradU&eX)*eX
                           + magSqr(gradU&eY)*eY
                           + magSqr(gradU&eZ)*eZ;
    
    // gradient of solids volume fraction
    volVectorField gradAlpha  = fvc::grad(alpha);
    
    // get turbulent kinetic energy of continuous-phase
    const volVectorField& kD(mesh_.lookupObject<volVectorField>
                                     ("k." + fluid.otherPhase(phase_).name()));
    
    const cellList& cells = mesh_.cells();
    
    // get drag coefficient
    volScalarField beta
    (
        fluid.lookupSubModel<dragModel>
        (
            fluid.otherPhase(phase_),
            phase_
        ).K()
    );
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
    
    if (dynamicAdjustment_) {
        // initialize
        volScalarField xiPhiGt1 = xiPhiG_&eX;
        volScalarField xiPhiGt2 = xiPhiG_&eY;
        volScalarField xiPhiGt3 = xiPhiG_&eZ;
        volScalarField U1 = U&eX;
        volScalarField U2 = U&eY;
        volScalarField U3 = U&eZ;
        // compute correlation between alpha.dispersed with U.continuous
        xiPhiGt1 = (
                      filter_(alpha*U1)
                    - filter_(alpha)*filter_(U1)
                   )/
                   (
                      sqrt(max(mag(filter_(sqr(alpha)))-sqr(filter_(alpha)),sqr(residualAlpha_)))*
                      sqrt(max(
                          mag(filter_(sqr(U1))
                        - sqr(filter_(U1))),uSmall)
                      )
                   );
        xiPhiGt2 = (
                      filter_(alpha*U2)
                    - filter_(alpha)*filter_(U2)
                   )/
                   (
                      sqrt(max(mag(filter_(sqr(alpha)))-sqr(filter_(alpha)),sqr(residualAlpha_)))*
                      sqrt(max(
                          mag(filter_(sqr(U2))
                        - sqr(filter_(U2))),uSmall)
                      )
                   );
        xiPhiGt3 = (
                      filter_(alpha*U3)
                    - filter_(alpha)*filter_(U3)
                   )/
                   (
                      sqrt(max(mag(filter_(sqr(alpha)))-sqr(filter_(alpha)),sqr(residualAlpha_)))*
                      sqrt(max(
                          mag(filter_(sqr(U3))
                        - sqr(filter_(U3))),uSmall)
                      )
                   );
        // smooth correlation coefficient
        xiPhiGt1.max(-0.99);
        xiPhiGt1.min(0.99);
        xiPhiGt2.max(-0.99);
        xiPhiGt2.min(0.99);
        xiPhiGt3.max(-0.99);
        xiPhiGt3.min(0.99);
        // negative sign since xiPhiG is computed from cont. phase volume fraction
        xiPhiG_ = - filter_(xiPhiGt1*eX + xiPhiGt2*eY + xiPhiGt3*eZ);

        // compute correlation coefficient between
        volScalarField magUd = mag(Ud_);
        volScalarField magUc = mag(U);
        
        volScalarField xiGSt = (filter_(magUc*magUd)-filter_(magUc)*filter_(magUd))/
                    (
                        sqrt(max(mag(filter_(sqr(magUc))-sqr(filter_(magUc))),sqr(uSmall)))
                      * sqrt(max(mag(filter_(sqr(magUd))-sqr(filter_(magUd))),sqr(uSmall)))
                     );
        // smooth correlation coefficient
        xiGS_ = filter_(xiGSt);
        xiGS_.max(-0.99);
        xiGS_.min(0.99);
        
        // Currently no dynamic procedure for Cmu and Ceps
        // Set Cmu
        Cmu_ = CmuScalar_;
        // Set Ceps
        Ceps_   = CepsScalar_;
    } else {
        // the sign of xiPhiG should be opposite to the slip velocity
        xiPhiG_ = xiPhiContScalar_*(sign(uSlip&eX)*eX + sign(uSlip&eY)*eY + sign(uSlip&eY)*eY);
        xiGS_   = xiGSScalar_;
        Cmu_    = CmuScalar_;
        Ceps_   = CepsScalar_;
    }

    // compute grid size
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
    
    // compute mixing length
    volScalarField lm = Cmu_*deltaF_;
    lm.max(lSmall.value());

    // Compute k_
    if (!equilibrium_) {
        fv::options& fvOptions(fv::options::New(mesh_));

        // Construct the transport equation for k
        // --> Stefanie
        fvVectorMatrix kEqn
        (
            fvm::ddt(alpha, rho, k_)
          + fvm::div(alphaRhoPhi, k_)
          - fvc::Sp(fvc::ddt(alpha, rho) + fvc::div(alphaRhoPhi), k_)
         ==
            fvOptions(alpha, rho, k_)
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
        // Schneiderbauer (2017), equ. (56)
        volVectorField kold = k_+ kSmall*eSum ;
        forAll(cells,cellI)
        {
            for (int i=0; i<3; i++) {
                k_[cellI].component(i) =
                    sqr(
                         - betaA[cellI]*lm[cellI]
                         + Foam::sqrt(
                              sqr(betaA[cellI]*lm[cellI])
                            + 2.0 * lm[cellI]
                            * Foam::max(
                                 lm[cellI]*SijSij[cellI].component(i)
                               + betaA[cellI] * xiGS_[cellI]*Foam::sqrt(kD[cellI].component(i))
                               - KdUdrift[cellI].component(i)*uSlip[cellI].component(i)
                                  /(alpha[cellI]*rho[cellI]*Foam::sqrt(kold[cellI].component(i)))
                                , 0.
                              )
                           )
                        ) / sqr(Ceps_[cellI]);
            }
        }
        // limit k_
        k_ *= min(mag(k_),maxK_)/max(mag(k_),sqr(uSmall));
        
        k_.correctBoundaryConditions();
    }
    
    //- compute variance of solids volume fraction
    volScalarField denom = fvc::div(U) + Cmu_ * Ceps_ * sqrt(k_ & eSum)/lm;
    volScalarField signDenom = sign(denom);
    denom.max(kSmall.value());
    
    alphaP2Mean_ =   4.0 * (xiPhiG_ & xiPhiG_) *
                     (
                          (k_ & eX)*sqr(gradAlpha & eX)
                        + (k_ & eY)*sqr(gradAlpha & eY)
                        + (k_ & eZ)*sqr(gradAlpha & eZ)
                      ) / sqr(denom) *  signDenom;
    alphaP2Mean_.max(0);
    alphaP2Mean_ = min(alphaP2Mean_, alpha*(1.0 - alpha));
    // compute nut_ (Schneiderbauer, 2017; equ. (34))
    nut_ = pos((scalar(1.0) - alpha) - residualAlpha_)*alpha*sqrt(k_ & eSum)*lm;
    
    // Limit viscosity and add frictional viscosity
    nut_.min(maxNut_);
    
    Info<< "SA-TFM (continuous Phase):" << nl
    << "    max(nut) = " << max(nut_).value() << nl
    << "    max(k_)  = " << max(k_ & eSum).value() << endl;
}


// ************************************************************************* //

#include "createRDeltaT.H"

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"


Info<< "Reading transportProperties\n" << endl;
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();


// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();


// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);


// Construct incompressible turbulence model
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, mixture)
);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

#include "createMRF.H"
#include "createFvOptions.H"

// fields required to close SGS terms
Info<< "Defining SGS extension\n" << endl;

volScalarField Cmu_
(
     IOobject
     (
          "Cmu",
          runTime.timeName(),
          mesh,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
     ),
     mesh,
     dimensionedScalar("value", dimensionSet(0,0,0,0,0,0,0), 0.3)
);

dimensionedScalar Csigma_
(
     "Csigma",
     dimensionSet(0,0,0,0,0,0,0),
     2.0
);

volScalarField nutSigma_
(
     IOobject
     (
          "nutSigma",
          runTime.timeName(),
          mesh,
          IOobject::MUST_READ,
          IOobject::AUTO_WRITE
     ),
     mesh
);

volScalarField Cst_
(
     IOobject
     (
          "Cst",
          runTime.timeName(),
          mesh,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
     ),
     mesh,
     dimensionedScalar("value", dimensionSet(0,0,0,0,0,0,0), 0)
);

volScalarField deltaF_
(
     IOobject
     (
          "deltaF",
          runTime.timeName(),
          mesh,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
     ),
     mesh,
     dimensionedScalar("1", dimensionSet(0,1,0,0,0,0,0), 1.0)
);

volScalarField cellVolume
(
    IOobject
    (
        "cellVolume",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("one", dimensionSet(0,3,0,0,0,0,0), 1)
);

Info<< "Computing grid spacing\n" << endl;

cellVolume.ref() = mesh.V();
deltaF_ = pow(cellVolume,1.0/3.0);

volScalarField wD = wallDist(mesh).y();
deltaF_ = neg(wD)*deltaF_ + pos(wD)*min(deltaF_,wD);


volScalarField corrSurfaceTensionForce_
(
    IOobject
    (
        "corrSurfaceTensionForce",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("one", dimensionSet(0,0,0,0,0,0,0), 1)
);

simpleTestFilter filter_(mesh);
simpleFilter filterS_(mesh);
laplaceFilter filterL_(mesh,24.0);

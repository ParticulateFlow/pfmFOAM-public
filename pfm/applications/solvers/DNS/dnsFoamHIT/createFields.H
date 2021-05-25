#include "readTransportProperties.H"

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
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

volVectorField term1g
(
    IOobject
    (
        "vortTerm1.air",
        U.time().timeName(),
        U.mesh(),
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimensionSet(0, 0, -2, 0, 0), vector(0,0,0))
);
volVectorField term2g
(
    IOobject
    (
        "vortTerm2.air",
        U.time().timeName(),
        U.mesh(),
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimensionSet(0, 0, -2, 0, 0), vector(0,0,0))
);
volVectorField term5g
(
    IOobject
    (
        "vortTerm5.air",
        U.time().timeName(),
        U.mesh(),
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimensionSet(0, 0, -2, 0, 0), vector(0,0,0))
);


#include "createPhi.H"

mesh.setFluxRequired(p.name());

#include "readTurbulenceProperties.H"

/*
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
*/

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

volVectorField deltaU
(
    IOobject
    (
        "deltaU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    U
);

volVectorField deltaUEvolved
(
    IOobject
    (
        "deltaUEvolved",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    U
);

volVectorField URef
(
    IOobject
    (
        "URef",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    U
);

volVectorField URefEvolved
(
    IOobject
    (
        "URefEvolved",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    U
);

IOdictionary dbProperties
(
    IOobject
    (
        "dataBaseProperties",
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

word dataBaseName(dbProperties.lookupOrDefault<word>("dataBaseName","dataBase"));
word URefStateName(dbProperties.lookupOrDefault<word>("UName","U"));
word exactSolutionFieldName(dbProperties.lookupOrDefault<word>("exactSolutionFieldName","none"));
bool compareToExactSolution = false;
if (exactSolutionFieldName != "none") compareToExactSolution = true;
autoPtr<volVectorField> UExactSolution(NULL);

#include "createFvOptions.H"

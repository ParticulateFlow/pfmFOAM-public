
Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0,2,-2,0,0,0,0), 0.0)
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

volScalarField divU
(
    IOobject
    (
        "divU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::div(U)
);
divU.write();

volVectorField convectiveTermLinear
(
    IOobject
    (
        "convectiveTermLinear",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("zero",dimensionSet(0,1,-2,0,0,0,0), vector::zero)
);

volVectorField convectiveTermQuadratic
(
    IOobject
    (
        "convectiveTermQuadratic",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("zero",dimensionSet(0,1,-2,0,0,0,0), vector::zero)
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
bool correctP(dbProperties.lookupOrDefault<bool>("correctP",false));

if (correctP && !p.headerOk())
{
    FatalError << "Provide pressure field to apply pressure correction." << abort(FatalError);
}

scalar minKuu(dbProperties.lookupOrDefault<scalar>("minKuu",1e-5));

label pRefCell = 0;
scalar pRefValue = 0.0;

#include "createFvOptions.H"

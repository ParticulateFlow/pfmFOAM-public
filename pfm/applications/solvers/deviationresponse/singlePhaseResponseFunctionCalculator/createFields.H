Info<< "Reading field X_pu\n" << endl;
volScalarField X_pu
(
    IOobject
    (
        "X_pu",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field X_uu\n" << endl;
volVectorField X_uu
(
    IOobject
    (
        "X_uu",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volVectorField divX_uu
(
    IOobject
    (
        "divX_uu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimensionSet(0,-4,0,0,0,0,0), vector::zero)
);

surfaceScalarField phiX
(
    IOobject
    (
        "phiX",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fvc::flux(X_uu)
);

volScalarField RHS
(
    IOobject
    (
        "RHS",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0,-4,-1,0,0,0,0), 0.0)
);

IOdictionary dbProperties
(
    IOobject
    (
        "dbProperties",
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

label nXuuCorr(dbProperties.lookupOrDefault<label>("nXuuCorr",1));
label minSourceElement(dbProperties.lookupOrDefault<label>("minSourceElement",0));
label maxSourceElement(dbProperties.lookupOrDefault<label>("maxSourceElement",1e9));
wordList patches(dbProperties.lookupOrDefault<wordList>("patches",wordList(0)));
labelList components(dbProperties.lookup("components"));
label numCmpt = components.size();
PtrList<volVectorField> X_uu_allCmpt(numCmpt);
for (int i = 0; i < numCmpt; i++)
{
    X_uu_allCmpt.set
    (
        i,
        new volVectorField(X_uu)
    );   
}
//List<Field<vector>> X_uu_allCmpt(numCmpt);
bool internalSource(dbProperties.lookupOrDefault<bool>("internalSource",true));
scalar minXuu(dbProperties.lookupOrDefault<scalar>("minXuu",1e-5));
scalar minXuu2 = minXuu*minXuu;

word dataBaseName(dbProperties.lookupOrDefault<word>("dataBase", "dataBase"));
word UdbName(dbProperties.lookupOrDefault<word>("U_db", "U"));
scalar tstart(dbProperties.lookupOrDefault<scalar>("tstart", 0.0));
scalar tolerance(dbProperties.lookupOrDefault<scalar>("tolerance", 1e-9));
Foam::Time dbTime(fileName(dataBaseName), "", "../system", "../constant", false);
instantList timeDirs(dbTime.times());
if (timeDirs.size() == 0)
{
    FatalError << "database " << dataBaseName << " does not exist or is empty\n" << abort(FatalError);
}

scalar deltaTdb = timeDirs[1].value() - timeDirs[0].value();
for (label i = 1; i < timeDirs.size() - 1; i++)
{
    if (Foam::mag(timeDirs[i+1].value() - timeDirs[i].value() - deltaTdb) > tolerance)
    {
        FatalError << "varying time step detected in database\n" << abort(FatalError);
    }
}

if (Foam::mag(deltaTdb - runTime.deltaTValue()) > tolerance)
{
    FatalError << "different time steps in database and solver\n" << abort(FatalError);
}

instantList::iterator it=timeDirs.begin();
while(true)
{
    if (Foam::mag(tstart - it->value()) < tolerance) break;
    if (it == timeDirs.end())
    {
        FatalError << "could not find start time " << tstart << " in database\n" << abort(FatalError);
    }
    it++;
}
instantList::iterator itStart = it;
dbTime.setTime(*itStart, itStart->value());

scalar tsim = runTime.endTime().value() - runTime.startTime().value();
if (timeDirs[timeDirs.size()-1].value() - it->value() < tsim - tolerance)
{
    FatalError << "database too short for intended simulation time\n" << abort(FatalError);
}

Info << "Reading field Udb\n" << endl;

volVectorField U_db
(
    IOobject
    (
        UdbName,
        dbTime.timePath(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::AUTO_WRITE
    ),
    mesh
);

surfaceScalarField phi_db
(
    IOobject
    (
        "phi_db",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::flux(U_db)
);


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(X_pu, piso.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(X_pu.name());


singlePhaseTransportModel laminarTransport(U_db, phi_db);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U_db, phi_db, laminarTransport)
);

#include "createFvOptions.H"

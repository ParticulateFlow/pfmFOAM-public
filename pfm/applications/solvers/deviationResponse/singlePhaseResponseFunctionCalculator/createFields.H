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

// testing
volVectorField X_uu1
(
    IOobject
    (
        "X_uu1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimensionSet(0,-3,0,0,0,0,0), vector::zero)
);

volVectorField X_uu2
(
    IOobject
    (
        "X_uu2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimensionSet(0,-3,0,0,0,0,0), vector::zero)
);

volScalarField X_pu1
(
    IOobject
    (
        "X_pu1",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0,-2,-1,0,0,0,0), 0.0)
);

volScalarField X_pu2
(
    IOobject
    (
        "X_pu2",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0,-2,-1,0,0,0,0), 0.0)
);
// testing done


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

IOdictionary tsProperties
(
    IOobject
    (
        "tsProperties",
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

label nXuuCorr(tsProperties.lookupOrDefault<label>("nXuuCorr",1));
label minSourceElement(tsProperties.lookupOrDefault<label>("minSourceElement",0));
label maxSourceElement(tsProperties.lookupOrDefault<label>("maxSourceElement",1e9));
wordList patches(tsProperties.lookupOrDefault<wordList>("patches",wordList(0)));
labelList components(tsProperties.lookup("components"));
label numCmpt = components.size();

OFstream OS_settings("patches");
OS_settings << patches << endl;

PtrList<volVectorField> X_uu_allCmpt(numCmpt);
for (int i = 0; i < numCmpt; i++)
{
    X_uu_allCmpt.set
    (
        i,
        new volVectorField(X_uu)
    );
}

scalar minXuu(tsProperties.lookupOrDefault<scalar>("minXuu",1e-5));
scalar minXuu2 = minXuu*minXuu;


enum Modes {internal, boundary, integrated};
word selectedMode(tsProperties.lookup("mode"));
Modes mode;
if (selectedMode == "internal") mode = internal;
else if (selectedMode == "boundary") mode = boundary;
else if (selectedMode == "integrated")
{
    mode = integrated;
    // instead of looping over all origin cells, do calculation just once
    minSourceElement = 0;
    maxSourceElement = 1;
    // keep all cells regardless of norm(Xuu)
    minXuu2 = -1.0;
}
else FatalError << "no known mode specified\n" << abort(FatalError);

word timeSeriesName(tsProperties.lookupOrDefault<word>("timeSeries", "timeSeries"));
word UtsName(tsProperties.lookupOrDefault<word>("U_ts", "U"));
scalar tstart(tsProperties.lookupOrDefault<scalar>("tstart", 0.0));
scalar tolerance(tsProperties.lookupOrDefault<scalar>("tolerance", 1e-9));
Foam::Time tsTime(fileName(timeSeriesName), "", "../system", "../constant", false);
instantList timeDirs(tsTime.times());
if (timeDirs.size() == 0)
{
    FatalError << "time series " << timeSeriesName << " does not exist or is empty\n" << abort(FatalError);
}

scalar deltaTts = timeDirs[1].value() - timeDirs[0].value();
for (label i = 1; i < timeDirs.size() - 1; i++)
{
    if (Foam::mag(timeDirs[i+1].value() - timeDirs[i].value() - deltaTts) > tolerance)
    {
        FatalError << "varying time step detected in database\n" << abort(FatalError);
    }
}

if (Foam::mag(deltaTts - runTime.deltaTValue()) > tolerance)
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
tsTime.setTime(*itStart, itStart->value());

scalar tsim = runTime.endTime().value() - runTime.startTime().value();
if (timeDirs[timeDirs.size()-1].value() - it->value() < tsim - tolerance)
{
    FatalError << "database too short for intended simulation time\n" << abort(FatalError);
}

Info << "Reading field Uts\n" << endl;

volVectorField U_ts
(
    IOobject
    (
        UtsName,
        tsTime.timePath(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::AUTO_WRITE
    ),
    mesh
);

surfaceScalarField phi_ts
(
    IOobject
    (
        "phi_ts",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::flux(U_ts)
);


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(X_pu, piso.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(X_pu.name());


singlePhaseTransportModel laminarTransport(U_ts, phi_ts);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U_ts, phi_ts, laminarTransport)
);

#include "createFvOptions.H"

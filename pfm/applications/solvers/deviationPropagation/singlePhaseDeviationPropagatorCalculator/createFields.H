Info<< "Reading field K_pu\n" << endl;
volScalarField K_pu
(
    IOobject
    (
        "K_pu",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field K_uu\n" << endl;
volVectorField K_uu
(
    IOobject
    (
        "K_uu",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField Receiver
(
    IOobject
    (
        "Receiver",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0,0,0,0,0,0,0), 0.0)
);

surfaceScalarField phiK
(
    IOobject
    (
        "phiK",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fvc::flux(K_uu)
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

label minSourceElement(tsProperties.lookupOrDefault<label>("minSourceElement",0));
label maxSourceElement(tsProperties.lookupOrDefault<label>("maxSourceElement",1e9));
wordList patches(tsProperties.lookupOrDefault<wordList>("patches",wordList(0)));
labelList components(tsProperties.lookup("components"));
label numCmpt = components.size();

OFstream OS_settings("patches");
OS_settings << patches << endl;

PtrList<volVectorField> K_uu_allCmpt(numCmpt);
for (int i = 0; i < numCmpt; i++)
{
    K_uu_allCmpt.set
    (
        i,
        new volVectorField(K_uu)
    );
}

scalar minKuu(tsProperties.lookupOrDefault<scalar>("minKuu",1e-5));


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
    // keep all cells regardless of norm(Kuu)
    minKuu = -1.0;
}
else FatalError << "no known mode specified\n" << abort(FatalError);

fileName timeSeriesName(tsProperties.lookupOrDefault<fileName>("timeSeries", "timeSeries"));
word UtsName(tsProperties.lookupOrDefault<word>("U_ts", "U"));
scalar tstart(tsProperties.lookupOrDefault<scalar>("tstart", 0.0));
scalar tolerance(tsProperties.lookupOrDefault<scalar>("tolerance", 1e-9));
Foam::Time tsTime(timeSeriesName, "", "../system", "../constant", false);
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

double dTimeStepRatio = deltaTts/runTime.deltaTValue() + 0.5;
label timeStepRatio = (int)dTimeStepRatio;

if (Foam::mag(deltaTts/runTime.deltaTValue() - timeStepRatio) > tolerance)
{
    FatalError << "time step in database not multiple of solver time step\n" << abort(FatalError);
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
setRefCell(K_pu, piso.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(K_pu.name());


singlePhaseTransportModel laminarTransport(U_ts, phi_ts);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U_ts, phi_ts, laminarTransport)
);

#include "createFvOptions.H"

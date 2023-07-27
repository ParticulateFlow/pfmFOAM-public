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
volVectorField Utot
(
    IOobject
    (
        "Utot",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    Utot
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

fileName    timeSeriesName      (tsProperties.lookupOrDefault<fileName>("timeSeries", "timeSeries"));
word        UtsName             (tsProperties.lookupOrDefault<word>("U_ts", "U"));
labelList   refStates           (tsProperties.lookupOrDefault<labelList>("refStates", labelList(0)));
scalar      tstart              (tsProperties.lookupOrDefault<scalar>("tstart",0.0));
scalar      tstep               (tsProperties.lookupOrDefault<scalar>("tstep",1.0));
scalar      tref = tstart + tstep * refStates[0];
scalar      tolerance           (tsProperties.lookupOrDefault<scalar>("tolerance", 1e-9));
label       resetRefStateEvery  (tsProperties.lookupOrDefault<label>("resetRefStateEvery",1));
(void)resetRefStateEvery;

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

label timeStepRatio(deltaTts/runTime.deltaTValue() + 0.5);

if (Foam::mag(deltaTts/runTime.deltaTValue() - timeStepRatio) > tolerance)
{
    FatalError << "time step in database not multiple of solver time step\n" << abort(FatalError);
}

instantList::iterator it=timeDirs.begin();
while(true)
{
    if (Foam::mag(tref - it->value()) < tolerance) break;
    if (it == timeDirs.end())
    {
        FatalError << "could not find start time " << tref << " in database\n" << abort(FatalError);
    }
    it++;
}
instantList::iterator itStart = it;
tsTime.setTime(*itStart, itStart->value());

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

U == Utot - U_ts;
U.oldTime() == U;

#include "createPhi.H"

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, piso.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());


singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

#include "createMRF.H"
#include "createFvOptions.H"

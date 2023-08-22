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

Application
    pisoFoam

Description
    Transient solver for incompressible, turbulent flow, using the PISO
    algorithm.

    Sub-models include:
    - turbulence modelling, i.e. laminar, RAS or LES
    - run-time selectable MRF and finite volume options, e.g. explicit porosity

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pisoControl.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    turbulence->validate();


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
    label subStepCounter = 0;
    label nthRefStateUpdate = 0;
    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        if (subStepCounter % timeStepRatio == 0)
        {
            if (subStepCounter % resetRefStateEvery == 0)
            {
                tref = tstart+tstep*refStates[nthRefStateUpdate];
                it=timeDirs.begin();
                while(true)
                {
                    if (Foam::mag(tref - it->value()) < tolerance) break;
                    if (it == timeDirs.end())
                    {
                        FatalError << "could not find start time " << tref << " in database\n" << abort(FatalError);
                    }
                    it++;
                }
                nthRefStateUpdate++;
            }
            tsTime.setTime(*it, it->value());
            Info << "tst = " << it->value() << endl;
            volVectorField U_ts_curr
            (
                IOobject
                (
                    UtsName,
                    tsTime.timePath(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            );
            U_ts == U_ts_curr;
            phi_ts == fvc::flux(U_ts);
            if (subStepCounter % resetRefStateEvery == 0)
            {
                U == Utot - U_ts;
                U.oldTime() == U;
                phi == fvc::flux(U);
            }
            it++;
        }

        #include "CourantNo.H"

        // Pressure-velocity PISO corrector
        {
            #include "UEqnExpl.H"

            // --- PISO loop
            while (piso.correct())
            {
                #include "pEqn.H"
            }
        }

        laminarTransport.correct();
        turbulence->correct();

        Utot == U + U_ts;

        runTime.write();

        subStepCounter++;

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

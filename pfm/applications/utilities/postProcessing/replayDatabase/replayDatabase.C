/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    replayDatabase


Description
    Read time series and step through it with predefined step size

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "functionObjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    functionObjectList funcObs(runTime);

    scalar t;


    funcObs.start();

    for(int timeI = 0; timeI < timeDirs.size(); timeI += samplestep)
    {
        dbTime.setTime(timeDirs[timeI], timeI);
        t = dbTime.value();
        runTime.setTime(t,timeI);
  //      if(t < startTime) continue;
  //      if(t > endTime) continue;
        Info << "time = " << t << ", time index = " << timeI << endl;


        for (int i=0; i<distDataScalarFieldNames_.size(); i++)
        {
            distDataScalarFields_.append
            (
                new volScalarField
                (
                    IOobject
                    (
                        distDataScalarFieldNames_[i],
                        dbTime.timePath(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                     ),
                     mesh
                )
           );
        }

        for (int i=0; i<distDataVectorFieldNames_.size(); i++)
        {
            distDataVectorFields_.append
            (
                new volVectorField
                (
                    IOobject
                    (
                        distDataVectorFieldNames_[i],
                        dbTime.timePath(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                     ),
                     mesh
                )
            );
        }

        funcObs.execute();

        distDataScalarFields_.clear();
        distDataVectorFields_.clear();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }
    funcObs.end();
   

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //

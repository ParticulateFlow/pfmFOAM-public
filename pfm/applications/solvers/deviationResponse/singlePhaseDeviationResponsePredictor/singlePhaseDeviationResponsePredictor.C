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
    singlePhaseDeviationResponsePredictor

Description
    Predicts the evolution of a single-phase flow based on deviation reponse
    functions.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dataBase.H"
#include "fieldNorm.H"
#include "responseFunctions.H"
#include "referenceStates.H"
#include "fvOptions.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"



    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
  
    dataBase db(mesh);
    db.init();
    label URefStateListIndex = db.referenceS().findRefStateListIndex("volVectorField",URefStateName);
    scalar predictionTimeStep = db.predictionTimeStep();
    runTime.setDeltaT(predictionTimeStep);

    OFstream distanceFile("distances");
    if (compareToExactSolution)
    {
        distanceFile << "# time\tinitial distance\terror of evolved ref solution\terror of complete prediction" << endl;
    }

    label nearestRefState = -1;
    bool firstStep = true;
    scalar initialDistance = -1.0;


    while (runTime.loop())
    {
        Info << "Time = " << runTime.timeName() << nl << endl;
        if (!firstStep)
        {
            U.write();
            URefEvolved.write();
            deltaUEvolved.write();
            if (compareToExactSolution)
            {
                UExactSolution.reset
                (
                    new volVectorField
                        (
                            IOobject
                            (
                                exactSolutionFieldName,
                                runTime.timeName(),
                                mesh,
                                IOobject::MUST_READ,
                                IOobject::NO_WRITE
                            ),
                            mesh
                        )
                );
                
                scalar distance1 = db.fieldN().fieldsDistance(UExactSolution,URefEvolved);
                scalar distance2 = db.fieldN().fieldsDistance(UExactSolution,U);

                distanceFile << runTime.timeName() << "\t" << initialDistance << "\t" << distance1 << "\t" << distance2 << endl;
            }
        }
        
        nearestRefState = db.findNearestRefState(U,URefStateListIndex);
        URef == db.referenceS().exportVolVectorField(URefStateListIndex,nearestRefState);
        if (firstStep)
        {
            initialDistance = db.fieldN().fieldsDistance(U,URef);
            Info << "Performing calculation for initial distance " << initialDistance << endl;
        }
        deltaU == U - URef;
        
        #include "convoluteDeviation.H"

        URefEvolved == db.referenceS().exportVolVectorEvolvedField(URefStateListIndex,nearestRefState);
        U = URefEvolved + deltaUEvolved;
        U.correctBoundaryConditions(); 

        URef.write();
        deltaU.write();

        firstStep = false;

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;


    }
    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //

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
    singlePhaseDeviationPropagationPredictor

Description
    Predicts the evolution of a single-phase flow based on deviation propagators.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dataBase.H"
#include "fieldNorm.H"
#include "deviationPropagators.H"
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
    label URefStateListIndex = db.exportReferenceStates().findRefStateListIndex("volVectorField",URefStateName);
    scalar predictionTimeStep = db.predictionTimeStep();
    runTime.setDeltaT(predictionTimeStep);

    OFstream distanceFile("distances");
    OFstream refStatesFile("refStates");
    if (compareToExactSolution)
    {
        distanceFile << "# time\tinitial distance\terror of evolved ref solution\terror of complete prediction" << endl;
    }

    label nearestRefState = -1;
    scalar initialDistance = -1;
    label prevTimeIndex = runTime.timeIndex();
    label presentTimeIndex = runTime.timeIndex();
    instant prevTimeInstant(runTime.value(),runTime.timeName());
    instant presentTimeInstant(runTime.value(),runTime.timeName());

    while (runTime.loop())
    {
        Info << "Time = " << runTime.timeName() << nl << endl;
        presentTimeIndex = runTime.timeIndex();
        presentTimeInstant.value() = runTime.value();
        presentTimeInstant.name() = runTime.timeName();

        nearestRefState = db.findNearestRefState(U,URefStateListIndex,initialDistance);
        refStatesFile << nearestRefState << " ";
        URef == db.exportReferenceStates().exportVolVectorField(URefStateListIndex,nearestRefState);
        deltaU == U - URef;

        convectiveTermQuadratic = fvc::div(linearInterpolate(deltaU) & mesh.Sf(),deltaU);

        if (runTime.writeTime())
        {
            runTime.setTime(prevTimeInstant,prevTimeIndex);
            URef.write();
            deltaU.write();
            convectiveTermLinear.write();
            convectiveTermQuadratic.write();
            runTime.setTime(presentTimeInstant,presentTimeIndex);
        }
        prevTimeIndex = presentTimeIndex;
        prevTimeInstant.value() = runTime.value();
        prevTimeInstant.name() = runTime.timeName();

        #include "convoluteDeviation.H"

        URefEvolved == db.exportReferenceStates().exportVolVectorEvolvedField(URefStateListIndex,nearestRefState);
        U = URefEvolved + deltaUEvolved;
        U.correctBoundaryConditions();
        if (correctP)
        {
            #include "correctP.H"
        }
        divU = fvc::div(U);

        runTime.write();

        distanceFile << runTime.timeName() << "\t" << initialDistance;
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

            volScalarField magUExact(mag(UExactSolution()));
            scalar normalization = fvc::domainIntegrate(magUExact).value();

            scalar distance1 = db.fieldN().fieldsDistance(UExactSolution(),URefEvolved,normalization);
            scalar distance2 = db.fieldN().fieldsDistance(UExactSolution(),U,normalization);

            distanceFile << "\t" << distance1 << "\t" << distance2;
        }
        distanceFile << endl;

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;


    }
    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //

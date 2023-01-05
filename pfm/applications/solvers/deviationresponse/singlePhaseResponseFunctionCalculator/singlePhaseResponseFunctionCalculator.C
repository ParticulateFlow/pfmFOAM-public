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
    singlePhaseResponseFunctionCalculator

Description
    Calculates the deviation response function for incompressible single-phase
    flow.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pisoControl.H"
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
    #include "createElementList.H"
    #include "initContinuityErrs.H"

 //   turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    fileName nameC, nameX;
    if (internalSource)
    {
        nameC = "targetCells_internal";
        nameX = "X_uu_internal";
    }
    else
    {
        nameC = "targetCells_boundary";
        nameX = "X_uu_boundary";
    }
    OFstream OS_C(nameC);
    OFstream OS_X(nameX);

    maxSourceElement = (sourceElements.size() > maxSourceElement) ? maxSourceElement : sourceElements.size();
    dimensionedVector zeroVec("zero",dimensionSet(0,-3,0,0,0),vector::zero);
    
    const volVectorField::Boundary& Xbf = X_uu.boundaryField();
    for (label sourceElement = minSourceElement; sourceElement < maxSourceElement; sourceElement++)
    {
        forAll(components,c)
        {
            label cmpt = components[c];
            // reset simulation and database time
            runTime.setTime(runTime.startTime(),runTime.startTimeIndex());
            it = itStart;

            X_uu == zeroVec;

            if (globalNumbering.isLocal(sourceElements[sourceElement]))
            {
                label sourceElementLocalID = globalNumbering.toLocal(sourceElements[sourceElement]);
                if (internalSource)
                {
                    Pout << "\nStarting time loop for internal source element " << sourceElements[sourceElement] << " with local ID " << sourceElementLocalID << endl;
                    X_uu[sourceElementLocalID].component(cmpt) = 1.0 / mesh.V()[sourceElementLocalID];
                }
                else
                {
                    label faceI = faceIDperPatch[sourceElementLocalID];
                    label patchI = patchOwningFace[sourceElementLocalID];
                    Pout << "\nStarting time loop for boundary source element " << sourceElements[sourceElement] << " with local ID " << faceI << " on patch " << boundaryMesh.names()[patchI] << endl;
                    // check BC on patch and decide what value needs to be set
                    if (isA<fixedValueFvPatchVectorField>(Xbf[patchI]))
                    {
                        X_uu.boundaryFieldRef()[patchI][faceI].component(cmpt) = 1.0 / mesh.boundary()[patchI].magSf()[faceI];
                    }
                    else if (isA<zeroGradientFvPatchVectorField>(Xbf[patchI]))
                    {
                        continue;
                    }
                    else
                    {
                        FatalError << "velocity boundary condition has to be fixed value or zero gradient, but found different condition on patch " << boundaryMesh.names()[patchI] << abort(FatalError);
                    }
                }

            }

            // Barrier ?

            scalar CoNum = 0.0;
            scalar deltaT = 0.0;

            while (runTime.loop())
            {
                Info<< "Time = " << runTime.timeName() << nl << endl;

                deltaT = X_uu.mesh().time().deltaTValue();
                dbTime.setTime(*it, it->value());
                Info << "dbt = " << it->value() << endl;
                volVectorField U_db_curr
                (
                    IOobject
                    (
                        UdbName,
                        dbTime.timePath(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh
                );
                U_db == U_db_curr;
                phi_db = fvc::flux(U_db);

                scalarField sumPhi
                (
                    fvc::surfaceSum(mag(phiX))().primitiveField()
                );
                CoNum = 0.5*gMax(sumPhi/mesh.V().field())*deltaT;
                Info<< "Courant Number max: " << CoNum << endl;

                // Pressure-velocity PISO corrector
                {
                    #include "XuuEqn.H"

                    // --- PISO loop
                    while (piso.correct())
                    {
                        #include "pEqn.H"
                    }
                }

                laminarTransport.correct();
                turbulence->correct();

                runTime.write();

                it++;

                Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                    << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                    << nl << endl;
            }
            X_uu_allCmpt[c] = X_uu;
        }
        // TODO: problem with evaluation of divX_uu
        forAll(components,c)
        {
            label cmpt = components[c];
            divX_uu.component(cmpt) = fvc::div(X_uu_allCmpt[c]);
        }
        divX_uu.write();
        #include "reconstructAndWrite.H"
    }
    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //

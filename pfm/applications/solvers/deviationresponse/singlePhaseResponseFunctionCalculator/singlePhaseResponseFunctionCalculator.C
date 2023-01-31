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
    if (mode == internal)
    {
        nameC = "targetCells_internal";
        nameX = "X_uu_internal";
    }
    else if (mode == boundary)
    {
        nameC = "targetCells_boundary";
        nameX = "X_uu_boundary";
    }
    else
    {
        nameC = "";
        nameX = "X_uu_integrated";
    }
    OFstream OS_C(nameC);
    OFstream OS_X(nameX);

    maxSourceElement = (sourceElements.size() > maxSourceElement) ? maxSourceElement : sourceElements.size();
    dimensionedVector defaultVec("zero",dimensionSet(0,-3,0,0,0),vector::zero);
    
    const volVectorField::Boundary& Xbf = X_uu.boundaryField();
    for (label sourceElement = minSourceElement; sourceElement < maxSourceElement; sourceElement++)
    {
        forAll(components,c)
        {
            label cmpt = components[c];
            // reset simulation and database time
            runTime.setTime(runTime.startTime(),runTime.startTimeIndex());
            it = itStart;
            
            if (mode == integrated)
            {
                defaultVec.value() = vector::zero;
                defaultVec.value().component(cmpt) = 1.0;
            }

            X_uu == defaultVec;

            if (mode != integrated)
            {
                if (globalNumbering.isLocal(sourceElements[sourceElement]))
                {
                    label sourceElementLocalID = globalNumbering.toLocal(sourceElements[sourceElement]);
                    if (mode == internal)
                    {
                        Pout << "\nStarting time loop for internal source element " << sourceElements[sourceElement] << " with local ID " << sourceElementLocalID << endl;
                        Pout << "\n 1/volume = " << 1.0 / mesh.V()[sourceElementLocalID] << endl;
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
            }

            // Barrier ?

            scalar CoNum = 0.0;
            dimensionedScalar deltaT("zero",dimensionSet(0,0,1,0,0,0,0),0.0);

            while (runTime.loop())
            {
                Info<< "Time = " << runTime.timeName() << nl << endl;
// testing
                /*
label sourceElementLocalID = globalNumbering.toLocal(sourceElements[sourceElement]);
label faceI = faceIDperPatch[sourceElementLocalID];
label patchI = patchOwningFace[sourceElementLocalID];
X_uu.boundaryFieldRef()[patchI][faceI].component(cmpt) = runTime.timeOutputValue() / mesh.boundary()[patchI].magSf()[faceI];
*/
// testing over
                deltaT.value() = X_uu.mesh().time().deltaTValue();
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
                phi_ts = fvc::flux(U_ts);

                scalarField sumPhi
                (
                    fvc::surfaceSum(mag(phiX))().primitiveField()
                );
                CoNum = 0.5*gMax(sumPhi/mesh.V().field())*deltaT.value();
                Info<< "Courant Number max: " << CoNum << endl;

                // Pressure-velocity PISO corrector
                {
                    #include "XuuEqn.H"

                    // --- PISO loop
                    while (piso.correct())
                    {
                        #include "pEqn.H"
                        //#include "pTEqn.H"
                    }
                }

                laminarTransport.correct();
                turbulence->correct();
                
                // testing
                if (c==0)
                {
                    X_uu1 = X_uu;
                    X_uu1.write();
                    
                    X_pu1 = X_pu;
                    X_pu1.write();
                }
                else
                {
                    X_uu2 = X_uu;    
                    X_uu2.write();

                    X_pu2 = X_pu;
                    X_pu2.write();
                }
                // testing done

                runTime.write();

                it++;

                Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                    << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                    << nl << endl;
            }
            X_uu_allCmpt[c] = X_uu;
        }
        forAll(components,c)
        {
            label cmpt = components[c];
            volScalarField divXuuC(fvc::div(X_uu_allCmpt[c]));
            forAll(X_uu,cellI)
            {
                divX_uu[cellI].component(cmpt) = divXuuC[cellI];
            }
        }
        scalar norm1 = 0.0;
        forAll(X_uu, cellI)
        {
            norm1 += magSqr(divX_uu[cellI]);
        }
        Info << "abs(div(X)) = " << norm1 << endl;
        divX_uu.write();
        #include "reconstructAndWrite.H"
    }
    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //

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
    singlePhaseDeviationPropagatorCalculator

Description
    Calculates the deviation propagator for incompressible single-phase
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
    fileName nameC, nameK;
    if (mode == internal)
    {
        nameC = "targetCells_internal";
        nameK = "K_uu_bySenders_internal";
    }
    else if (mode == boundary)
    {
        nameC = "targetCells_boundary";
        nameK = "K_uu_bySenders_boundary";
    }
    else
    {
        nameC = "";
        nameK = "K_uu_bySenders_integrated";
    }
    OFstream OS_C(nameC);
    OFstream OS_K(nameK);

    maxSourceElement = (sourceElements.size() > maxSourceElement) ? maxSourceElement : sourceElements.size();
    dimensionedVector defaultVec("zero",dimensionSet(0,-3,0,0,0),vector::zero);

    scalar normalization = 1.0;

    volVectorField::Boundary& Kbf = K_uu.boundaryFieldRef();
    for (label sourceElement = minSourceElement; sourceElement < maxSourceElement; sourceElement++)
    {
        Receiver == 0.0;
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

            K_uu == defaultVec;
            K_uu_allCmpt[c] = K_uu;

            if (mode != integrated)
            {
                // make sure any previously set gradients on the boundary are reset to zero
                for (int i = 0; i < patches.size(); i++)
                {
                    label patchI = boundaryMesh.findPatchID(patches[i]);
                    if (isA<fixedGradientFvPatchVectorField>(Kbf[patchI]))
                    {
                        fixedGradientFvPatchVectorField& Kbf_fixedGradient = dynamic_cast<fixedGradientFvPatchVectorField &>(Kbf[patchI]);
                        Kbf_fixedGradient.gradient() = vector::zero;
                    }
                }

                if (globalNumbering.isLocal(sourceElements[sourceElement]))
                {
                    label sourceElementLocalID = globalNumbering.toLocal(sourceElements[sourceElement]);
                    if (mode == internal)
                    {
                        normalization = mesh.V()[sourceElementLocalID];
                        Pout << "\nStarting time loop for internal source element " << sourceElements[sourceElement] << " with local ID " << sourceElementLocalID << endl;
                        Pout << "\n 1/volume = " << 1.0 / normalization << endl;
                        K_uu[sourceElementLocalID].component(cmpt) = 1.0 / normalization;

                        // adjacent to fixed gradient boundary?
                        if (adjacentFaceID[sourceElementLocalID] >= 0)
                        {
                            label faceI = adjacentFaceID[sourceElementLocalID];
                            label patchI = adjacentPatchID[sourceElementLocalID];
                       //     scalar normalization2 = mesh.boundary()[patchI].magSf()[faceI];
                            fixedGradientFvPatchVectorField& Kbf_fixedGradient = dynamic_cast<fixedGradientFvPatchVectorField &>(Kbf[patchI]);
                            Kbf_fixedGradient.gradient()[faceI].component(cmpt) = -1.0 / normalization * Kbf[patchI].patch().deltaCoeffs()[faceI];
                        }
                    }
                    else
                    {
                        label faceI = faceIDperPatch[sourceElementLocalID];
                        label patchI = patchOwningFace[sourceElementLocalID];
                        normalization = mesh.boundary()[patchI].magSf()[faceI];
                        Pout << "\nStarting time loop for boundary source element " << sourceElements[sourceElement] << " with local ID " << faceI << " on patch " << boundaryMesh.names()[patchI] << endl;
                        Pout << "\n 1/area = " << 1.0 / normalization << endl;
                        // check BC on patch and decide what value needs to be set
                        if (isA<fixedValueFvPatchVectorField>(Kbf[patchI]))
                        {
                            Kbf[patchI][faceI].component(cmpt) = 1.0 / normalization;
                        }
                        else if (isA<fixedGradientFvPatchVectorField>(Kbf[patchI]))
                        {
                            fixedGradientFvPatchVectorField& Kbf_fixedGradient = dynamic_cast<fixedGradientFvPatchVectorField &>(Kbf[patchI]);
                            Kbf_fixedGradient.gradient()[faceI].component(cmpt) = 1.0 / normalization * Kbf[patchI].patch().deltaCoeffs()[faceI];
                        }
                        else
                        {
                            FatalError << "velocity boundary condition has to be fixed value or fixed gradient, but found different condition on patch " << boundaryMesh.names()[patchI] << abort(FatalError);
                        }
                    }
                }
            }

            // Barrier ?

            scalar CoNum = 0.0;
            dimensionedScalar deltaT("zero",dimensionSet(0,0,1,0,0,0,0),0.0);
            deltaT.value() = K_uu.mesh().time().deltaTValue();
            label subStepCounter = 0;

            while (runTime.loop())
            {
                Info<< "Time = " << runTime.timeName() << nl << endl;
                if (subStepCounter % timeStepRatio == 0)
                {
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
                    it++;
                }
                scalarField sumPhi
                (
                    fvc::surfaceSum(mag(phiK))().primitiveField()
                );
                CoNum = 0.5*gMax(sumPhi/mesh.V().field())*deltaT.value();
                Info<< "Courant Number max: " << CoNum << endl;

                // Pressure-velocity PISO corrector
                {
                    #include "KuuEqnExpl.H"
                    #include "KpuEqnExpl.H"
                }

                laminarTransport.correct();
                turbulence->correct();

                runTime.write();

                subStepCounter++;

                Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                    << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                    << nl << endl;
            }
            K_uu_allCmpt[c] = K_uu;
        }
        #include "reconstructAndWrite.H"
    }

    tsTime.setTime(*it, it->value());
    OFstream endTime("timeEvolvedRefState");
    endTime << it->value() << endl;
    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //

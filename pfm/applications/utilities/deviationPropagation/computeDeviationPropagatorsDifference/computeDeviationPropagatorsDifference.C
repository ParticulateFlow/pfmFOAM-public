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
    computeDeviationPropagatorsDifference

Description
    Compute difference between pairs of deviation propagators.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "OFstream.H"
#include "dataBase.H"
#include "fieldNorm.H"
#include "deviationPropagators.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    if (!Pstream::master())
    {
        FatalError <<"current implementation for serial runs only\n" << abort(FatalError);
    }

    dataBase db1(mesh,"dataBase1");
    dataBase db2(mesh,"dataBase2");
    db1.init();
    db2.init();
    OFstream distanceFile("KuuDistances");
    distanceFile << "# refState distance" << endl;

    label numRefStates1 = db1.numRefStates();
    label numRefStates2 = db2.numRefStates();
    if (numRefStates1 != numRefStates2)
    {
        FatalError <<"different number of reference states in databases 1 and 2\n" << abort(FatalError);
    }
    scalar domainVol = gSum(mesh.V());

    for (int refState = 0; refState < numRefStates1; refState++)
    {
        word fieldName = "deltaK_uu_integrated_"+name(refState);
        volScalarField deltaK_uu_integrated
        (
            IOobject
            (
                fieldName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero",dimensionSet(0,0,0,0,0,0,0),0)
        );

        forAll(deltaK_uu_integrated, cellI)
        {
            deltaK_uu == dimensionedTensor("zero",dimensionSet(0,0,-3,0,0,0,0),tensor::zero);

            labelList &senderCells1 = db1.exportDeviationPropagators().senderCellIDs(refState,cellI);
            tensorList &K_uu1 = db1.exportDeviationPropagators().Kuu_internal(refState,cellI);

            labelList &senderCells2 = db2.exportDeviationPropagators().senderCellIDs(refState,cellI);
            tensorList &K_uu2 = db2.exportDeviationPropagators().Kuu_internal(refState,cellI);

            scalar norm1 = 0.0;
            scalar norm2 = 0.0;

            forAll(K_uu1, sender)
            {
                deltaK_uu[senderCells1[sender]] += K_uu1[sender];
                norm1 += Foam::sqrt(magSqr(K_uu1[sender])) * mesh.V()[senderCells1[sender]];
            }

            forAll(K_uu2, sender)
            {
                deltaK_uu[senderCells2[sender]] -= K_uu2[sender];
                norm2 += Foam::sqrt(magSqr(K_uu2[sender])) * mesh.V()[senderCells2[sender]];
            }

            labelList &senderBoundaryFaces1 = db1.exportDeviationPropagators().senderBoundaryFaceIDs(refState,cellI);
            tensorList &K_uu_boundary1 = db1.exportDeviationPropagators().Kuu_boundary(refState,cellI);
            labelList &faceIDperPatch1 = db1.faceIDperPatch();
            labelList &patchOwningFace1 = db1.patchOwningFace();

            labelList &senderBoundaryFaces2 = db2.exportDeviationPropagators().senderBoundaryFaceIDs(refState,cellI);
            tensorList &K_uu_boundary2 = db2.exportDeviationPropagators().Kuu_boundary(refState,cellI);
            labelList &faceIDperPatch2 = db2.faceIDperPatch();
            labelList &patchOwningFace2 = db2.patchOwningFace();

            label patchID;
            label faceID;
            forAll(K_uu_boundary1, sender)
            {
                patchID = patchOwningFace1[senderBoundaryFaces1[sender]];
                faceID = faceIDperPatch1[senderBoundaryFaces1[sender]];
                deltaK_uu.boundaryFieldRef()[patchID][faceID] += K_uu_boundary1[sender];
                norm1 += Foam::sqrt(magSqr(K_uu_boundary1[sender])) * mesh.magSf().boundaryField()[patchID][faceID];
            }

            forAll(K_uu_boundary2, sender)
            {
                patchID = patchOwningFace2[senderBoundaryFaces2[sender]];
                faceID = faceIDperPatch2[senderBoundaryFaces2[sender]];
                deltaK_uu.boundaryFieldRef()[patchID][faceID] -= K_uu_boundary2[sender];
                norm2 += Foam::sqrt(magSqr(K_uu_boundary2[sender])) * mesh.magSf().boundaryField()[patchID][faceID];
            }

            scalar distanceLocal = 0.0;
            scalar value = 0.0;
            forAll(deltaK_uu,cellJ)
            {
                value = magSqr(deltaK_uu[cellJ]);
                if (value >= VSMALL)
                {
                    distanceLocal += Foam::sqrt(value) * mesh.V()[cellJ];
                }
            }
            forAll(mesh.boundary(), patchJ) 
            {
                forAll (mesh.boundary()[patchJ],faceJ) 
                {
                    value = magSqr(deltaK_uu.boundaryField()[patchJ][faceJ]);
                    if (value >= VSMALL)
                    {
                        distanceLocal += Foam::sqrt(value) * mesh.magSf().boundaryField()[patchJ][faceJ];
                    }
                }
            }
            deltaK_uu_integrated[cellI] = distanceLocal / (Foam::sqrt(norm1) * Foam::sqrt(norm2));
        }
        scalar distance = fvc::domainIntegrate(deltaK_uu_integrated).value()/domainVol;
        distanceFile << refState << " " << distance << endl;
        deltaK_uu_integrated.write();
    }

    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //


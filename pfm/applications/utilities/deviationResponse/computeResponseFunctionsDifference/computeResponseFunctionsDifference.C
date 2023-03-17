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
    computeResponseFunctionsDifference

Description
    Compute difference between pairs of response functions.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "OFstream.H"
#include "dataBase.H"
#include "fieldNorm.H"
#include "responseFunctions.H"
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
    OFstream distanceFile("XuuDistances");
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
        word fieldName = "deltaX_uu_integrated_"+name(refState);
        volScalarField deltaX_uu_integrated
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

        forAll(deltaX_uu_integrated, cellI)
        {
            deltaX_uu == dimensionedTensor("zero",dimensionSet(0,0,-3,0,0,0,0),tensor::zero);

            labelList &senderCells1 = db1.responseF().senderCellIDs(refState,cellI);
            tensorList &X_uu1 = db1.responseF().Xuu_internal(refState,cellI);

            labelList &senderCells2 = db2.responseF().senderCellIDs(refState,cellI);
            tensorList &X_uu2 = db2.responseF().Xuu_internal(refState,cellI);

            forAll(X_uu1, sender)
            {
                deltaX_uu[senderCells1[sender]] += X_uu1[sender];
            }

            forAll(X_uu2, sender)
            {
                deltaX_uu[senderCells2[sender]] -= X_uu2[sender];
            }

            labelList &senderBoundaryFaces1 = db1.responseF().senderBoundaryFaceIDs(refState,cellI);
            tensorList &X_uu_boundary1 = db1.responseF().Xuu_boundary(refState,cellI);
            labelList &faceIDperPatch1 = db1.faceIDperPatch();
            labelList &patchOwningFace1 = db1.patchOwningFace();

            labelList &senderBoundaryFaces2 = db2.responseF().senderBoundaryFaceIDs(refState,cellI);
            tensorList &X_uu_boundary2 = db2.responseF().Xuu_boundary(refState,cellI);
            labelList &faceIDperPatch2 = db2.faceIDperPatch();
            labelList &patchOwningFace2 = db2.patchOwningFace();

            label patchID;
            label faceID;
            forAll(X_uu_boundary1, sender)
            {
                patchID = patchOwningFace1[senderBoundaryFaces1[sender]];
                faceID = faceIDperPatch1[senderBoundaryFaces1[sender]];
                deltaX_uu.boundaryFieldRef()[patchID][faceID] += X_uu_boundary1[sender];
            }

            forAll(X_uu_boundary2, sender)
            {
                patchID = patchOwningFace2[senderBoundaryFaces2[sender]];
                faceID = faceIDperPatch2[senderBoundaryFaces2[sender]];
                deltaX_uu.boundaryFieldRef()[patchID][faceID] -= X_uu_boundary2[sender];
            }

            scalar distanceLocal = 0.0;
            scalar value = 0.0;
            forAll(deltaX_uu,cellJ)
            {
                value = magSqr(deltaX_uu[cellJ]);
                if (value >= VSMALL)
                {
                    distanceLocal += Foam::sqrt(value) * mesh.V()[cellJ];
                }
            }
            forAll(mesh.boundary(), patchJ) 
            {
                forAll (mesh.boundary()[patchJ],faceJ) 
                {
                    value = magSqr(deltaX_uu.boundaryField()[patchJ][faceJ]);
                    if (value >= VSMALL)
                    {
                        distanceLocal += Foam::sqrt(value) * mesh.magSf().boundaryField()[patchJ][faceJ];
                    }
                }
            }

            deltaX_uu_integrated[cellI] = distanceLocal;
        }
        scalar distance = fvc::domainIntegrate(deltaX_uu_integrated).value()/domainVol;
        distanceFile << refState << " " << distance << endl;
        deltaX_uu_integrated.write();
    }

    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //


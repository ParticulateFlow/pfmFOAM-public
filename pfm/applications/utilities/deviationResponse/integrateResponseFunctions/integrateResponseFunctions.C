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
    integrateResponseFunctions

Description
    Integrates response functions over their second spatial argument.

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

    dataBase db(mesh);
    db.responseF().readIntegratedResponseFunctions();
    db.init();
    OFstream distanceFile("integratedXuuDistances");
    distanceFile << "# refState distance" << endl;

    label numRefStates = db.numRefStates();

    for (int refState = 0; refState < numRefStates; refState++)
    {
        runTime.setTime(refState,refState);

        forAll(X_uu_integrated, cellI)
        {
            X_uu_integrated_target[cellI] = db.responseF().Xuu_integrated(refState,cellI);

            labelList &senderCells = db.responseF().senderCellIDs(refState,cellI);
            tensorList &X_uu = db.responseF().Xuu_internal(refState,cellI);

            X_uu_integrated[cellI] = tensor::zero;
            scalar norm = 0.0;
            forAll(X_uu, sender)
            {
                norm = X_uu[sender].xx()*X_uu[sender].xx() +
                    X_uu[sender].xy()*X_uu[sender].xy() +
                    X_uu[sender].xz()*X_uu[sender].xz() +
                    X_uu[sender].yx()*X_uu[sender].yx() +
                    X_uu[sender].yy()*X_uu[sender].yy() +
                    X_uu[sender].yz()*X_uu[sender].yz() +
                    X_uu[sender].zx()*X_uu[sender].zx() +
                    X_uu[sender].zy()*X_uu[sender].zy() +
                    X_uu[sender].zz()*X_uu[sender].zz();

                if (Foam::sqrt(norm) * mesh.V()[senderCells[sender]] > minXuu)
                {
                    X_uu_integrated[cellI] += X_uu[sender] * mesh.V()[senderCells[sender]];
                }
            }

            labelList &senderBoundaryFaces = db.responseF().senderBoundaryFaceIDs(refState,cellI);
            tensorList &X_uu_boundary = db.responseF().Xuu_boundary(refState,cellI);
            labelList &faceIDperPatch = db.faceIDperPatch();
            labelList &patchOwningFace = db.patchOwningFace();

            label patchID;
            label faceID;
            norm = 0.0;
            forAll(X_uu_boundary, sender)
            {
                norm = X_uu_boundary[sender].xx()*X_uu_boundary[sender].xx() +
                    X_uu_boundary[sender].xy()*X_uu_boundary[sender].xy() +
                    X_uu_boundary[sender].xz()*X_uu_boundary[sender].xz() +
                    X_uu_boundary[sender].yx()*X_uu_boundary[sender].yx() +
                    X_uu_boundary[sender].yy()*X_uu_boundary[sender].yy() +
                    X_uu_boundary[sender].yz()*X_uu_boundary[sender].yz() +
                    X_uu_boundary[sender].zx()*X_uu_boundary[sender].zx() +
                    X_uu_boundary[sender].zy()*X_uu_boundary[sender].zy() +
                    X_uu_boundary[sender].zz()*X_uu_boundary[sender].zz();

                patchID = patchOwningFace[senderBoundaryFaces[sender]];
                faceID = faceIDperPatch[senderBoundaryFaces[sender]];
                if (Foam::sqrt(norm) * mesh.magSf().boundaryField()[patchID][faceID] > minXuu)
                {
                    X_uu_integrated[cellI] += X_uu_boundary[sender] * mesh.magSf().boundaryField()[patchID][faceID];
                }
            }
        }

        scalar distance = db.fieldN().fieldsDistance(X_uu_integrated,X_uu_integrated_target);
        distanceFile << refState << " " << distance << endl;
        X_uu_integrated.write();
        X_uu_integrated_target.write();
    }

    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //


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
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    dataBase db(mesh);
    db.init();
    
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
            forAll(X_uu, sender)
            {
                X_uu_integrated[cellI] += X_uu[sender] * mesh.V()[senderCells[sender]];
            }

            labelList &senderBoundaryFaces = db.responseF().senderBoundaryFaceIDs(refState,cellI);
            tensorList &X_uu_boundary = db.responseF().Xuu_boundary(refState,cellI);
            labelList &faceIDperPatch = db.faceIDperPatch();
            labelList &patchOwningFace = db.patchOwningFace();

            label patchID;
            label faceID;
            forAll(X_uu_boundary, sender)
            {
                patchID = patchOwningFace[sender];
                faceID = faceIDperPatch[sender];
                X_uu_integrated[cellI] += X_uu[sender] * mesh.magSf().boundaryField()[patchID][faceID];
            }
        }
        
        X_uu_integrated.write();
        X_uu_integrated_target.write();
    }

    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //


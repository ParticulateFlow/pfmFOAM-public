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
    integrateDeviationPropagators

Description
    Integrates deviation propagators over their second spatial argument.

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

    dataBase db(mesh);
    db.exportDeviationPropagators().readIntegratedDeviationPropagators();
    db.init();
    OFstream distanceFile("integratedKuuDistances");
    distanceFile << "# refState distance" << endl;

    label numRefStates = db.numRefStates();

    for (int refState = 0; refState < numRefStates; refState++)
    {
        runTime.setTime(refState,refState);

        scalar K_uu_integrated_target_sum = 0.0;
        forAll(K_uu_integrated, cellI)
        {
            K_uu_integrated_target[cellI] = db.exportDeviationPropagators().Kuu_integrated(refState,cellI);
            K_uu_integrated_target_sum += Foam::sqrt(magSqr(K_uu_integrated_target[cellI]))* mesh.V()[cellI];

            labelList &senderCells = db.exportDeviationPropagators().senderCellIDs(refState,cellI);
            tensorList &K_uu = db.exportDeviationPropagators().Kuu_internal(refState,cellI);

            K_uu_integrated[cellI] = tensor::zero;
            scalar norm = 0.0;
            forAll(K_uu, sender)
            {
                norm = K_uu[sender].xx()*K_uu[sender].xx() +
                    K_uu[sender].xy()*K_uu[sender].xy() +
                    K_uu[sender].xz()*K_uu[sender].xz() +
                    K_uu[sender].yx()*K_uu[sender].yx() +
                    K_uu[sender].yy()*K_uu[sender].yy() +
                    K_uu[sender].yz()*K_uu[sender].yz() +
                    K_uu[sender].zx()*K_uu[sender].zx() +
                    K_uu[sender].zy()*K_uu[sender].zy() +
                    K_uu[sender].zz()*K_uu[sender].zz();

                if (Foam::sqrt(norm) * mesh.V()[senderCells[sender]] > minKuu)
                {
                    K_uu_integrated[cellI] += K_uu[sender] * mesh.V()[senderCells[sender]];
                }
            }

            labelList &senderBoundaryFaces = db.exportDeviationPropagators().senderBoundaryFaceIDs(refState,cellI);
            tensorList &K_uu_boundary = db.exportDeviationPropagators().Kuu_boundary(refState,cellI);
            labelList &faceIDperPatch = db.faceIDperPatch();
            labelList &patchOwningFace = db.patchOwningFace();

            label patchID;
            label faceID;
            norm = 0.0;
            forAll(K_uu_boundary, sender)
            {
                norm = K_uu_boundary[sender].xx()*K_uu_boundary[sender].xx() +
                    K_uu_boundary[sender].xy()*K_uu_boundary[sender].xy() +
                    K_uu_boundary[sender].xz()*K_uu_boundary[sender].xz() +
                    K_uu_boundary[sender].yx()*K_uu_boundary[sender].yx() +
                    K_uu_boundary[sender].yy()*K_uu_boundary[sender].yy() +
                    K_uu_boundary[sender].yz()*K_uu_boundary[sender].yz() +
                    K_uu_boundary[sender].zx()*K_uu_boundary[sender].zx() +
                    K_uu_boundary[sender].zy()*K_uu_boundary[sender].zy() +
                    K_uu_boundary[sender].zz()*K_uu_boundary[sender].zz();

                patchID = patchOwningFace[senderBoundaryFaces[sender]];
                faceID = faceIDperPatch[senderBoundaryFaces[sender]];
                if (Foam::sqrt(norm) * mesh.magSf().boundaryField()[patchID][faceID] > minKuu)
                {
                    K_uu_integrated[cellI] += K_uu_boundary[sender] * mesh.magSf().boundaryField()[patchID][faceID];
                }
            }
        }

        reduce(K_uu_integrated_target_sum, sumOp<scalar>());
        scalar distance = db.fieldN().fieldsDistance(K_uu_integrated,K_uu_integrated_target,1.0);
        distance /= K_uu_integrated_target_sum;
        distanceFile << refState << " " << distance << endl;
        K_uu_integrated.write();
        K_uu_integrated_target.write();
    }

    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //


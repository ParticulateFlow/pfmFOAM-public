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
    illustrateDeviationPropagators

Description
    Writes deviation propagators for specific first argument.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "OFstream.H"
#include "dataBase.H"
#include "deviationPropagators.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    IOdictionary illustrateProperties
    (
        IOobject
        (
            "illustrateProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    label cellIMin(illustrateProperties.lookupOrDefault<label>("cellIMin",0));
    label cellIMax(illustrateProperties.lookupOrDefault<label>("cellIMax",0));
    label refStateMin(illustrateProperties.lookupOrDefault<label>("refStateMin",0));
    label refStateMax(illustrateProperties.lookupOrDefault<label>("refStateMax",0));

 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    if (!Pstream::master())
    {
        FatalError <<"current implementation for serial runs only\n" << abort(FatalError);
    }

    dataBase db(mesh);
    db.init();

    label numRefStates = db.numRefStates();

    for (label refState = refStateMin; refState <= refStateMax; refState++)
    {
        for (label cellI = cellIMin; cellI <= cellIMax; cellI++)
        {
            word fieldName = "K_uu_"+name(cellI)+"_"+name(refState);
            volTensorField K_uu_field
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
                dimensionedTensor("zero",dimensionSet(0,0,-3,0,0,0,0),tensor::zero)
            );

            labelList &senderCells = db.exportDeviationPropagators().senderCellIDs(refState,cellI);
            tensorList &K_uu = db.exportDeviationPropagators().Kuu_internal(refState,cellI);
            forAll(K_uu, sender)
            {
                K_uu_field[senderCells[sender]] = K_uu[sender];
            }

            labelList &senderBoundaryFaces = db.exportDeviationPropagators().senderBoundaryFaceIDs(refState,cellI);
            tensorList &K_uu_boundary = db.exportDeviationPropagators().Kuu_boundary(refState,cellI);
            labelList &faceIDperPatch = db.faceIDperPatch();
            labelList &patchOwningFace = db.patchOwningFace();

            label patchID;
            label faceID;
            forAll(K_uu_boundary, sender)
            {
                patchID = patchOwningFace[senderBoundaryFaces[sender]];
                faceID = faceIDperPatch[senderBoundaryFaces[sender]];
                K_uu_field.boundaryFieldRef()[patchID][faceID] = K_uu_boundary[sender];
            }

            K_uu_field.write();
        }
    }

    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //


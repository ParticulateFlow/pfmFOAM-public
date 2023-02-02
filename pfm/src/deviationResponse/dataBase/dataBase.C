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

\*---------------------------------------------------------------------------*/

#include "referenceStates.H"
#include "responseFunctions.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
dataBase::dataBase
(
    const fvMesh& mesh
)
:
    regIOobject
    (
        IOobject
        (
            "dataBase",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true
        )
    ),
    mesh_(mesh),
    dataBaseProperties_
    (
        IOobject
        (
            "dataBaseProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    dataBaseNames_(dataBaseProperties_.lookupOrDefault<wordList>("dataBases", wordList(1,"dataBase"))),
    numDataBases_(dataBaseNames_.size()),
    referenceStates_
    (
        referenceStates::New
        (
            dataBaseProperties_,
            *this
        )
    ),
    responseFunctions_
    (
        responseFunctions::New
        (
            dataBaseProperties_,
            *this
        )
    ),
    globalCellNumbering_(),
    globalBoundaryFaceNumbering_(),
    faceIDperPatch_(0),
    patchOwningFace_(0)
{
    globalCellNumbering_.set(new globalIndex(mesh.nCells()));
}




// * * * * * * * * * * * * * * * * Destructors  * * * * * * * * * * * * * * //
dataBase::~dataBase()
{}

// * * * * * * * * * * * * * * * public Member Functions  * * * * * * * * * * * * * //

const fvMesh& dataBase::mesh() const
{
    return mesh_;
}

void dataBase::init()
{
    wordList patches;
    IFstream IS("patches");
    IS >> patches;

    label locFaces = 0;
    for (int i = 0; i < patches.size(); i++)
    {
        label patchID = mesh_.boundaryMesh().findPatchID(patches[i]);
        locFaces += mesh_.boundaryMesh()[patchID].size();
    }
    faceIDperPatch_.resize(locFaces);
    patchOwningFace_.resize(locFaces);

    label fc = 0;
    for (int i = 0; i < patches.size(); i++)
    {
        label patchID = mesh_.boundaryMesh().findPatchID(patches[i]);

        forAll (mesh_.boundaryMesh()[patchID],faceI)
        {
            faceIDperPatch_[fc] = faceI;
            patchOwningFace_[fc] = patchID;
            fc++;
        }
    }

    globalBoundaryFaceNumbering_.set(new globalIndex(locFaces));

    responseFunctions_->readSenderIDs(dataBaseNames_);
    responseFunctions_->readResponseFunctions(dataBaseNames_);
}

referenceStates& dataBase::referenceS()
{
    return referenceStates_();
}

responseFunctions& dataBase::responseF()
{
    return responseFunctions_();
}

label dataBase::localFromGlobalCellID(label cellI)
{
    label id = -1;
    if (globalCellNumbering_().isLocal(cellI))
    {
        id = globalCellNumbering_().toLocal(cellI);
    }
    return id;
}

label dataBase::localFromGlobalBoundaryFaceID(label bFaceI)
{
    label id = -1;
    if (globalBoundaryFaceNumbering_().isLocal(bFaceI))
    {
        id = globalBoundaryFaceNumbering_().toLocal(bFaceI);
    }
    return id;
}



}

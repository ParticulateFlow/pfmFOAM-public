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

#include "fieldNorm.H"
#include "referenceStates.H"
#include "responseFunctions.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
dataBase::dataBase
(
    const fvMesh& mesh, word dbName
)
:
    regIOobject
    (
        IOobject
        (
            dbName,
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
            dbName+"Properties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    dataBaseNames_(dataBaseProperties_.lookupOrDefault<fileNameList>("dataBases", fileNameList(1,"dataBase"))),
    numDataBases_(dataBaseNames_.size()),
    numRefStates_(0),
    predictionTimeStep_(-1.0),
    fieldNorm_
    (
        fieldNorm::New
        (
            dataBaseProperties_,
            *this
        )
    ),
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
    // read reference states
    label numRefStates = -1;
    if (referenceStates_.valid())
    {
        numRefStates = referenceStates_->readReferenceStates(dataBaseNames_);
    }

    // in case only reference states are to be read, stop here
    if (dataBaseProperties_.lookup("responseFunctions") == "noResponseFunctions")
    {
        return;
    }
    
    // construct face addressing
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

    // read response functions
    numRefStates_ = responseFunctions_->readSenderIDs(dataBaseNames_);
    if (numRefStates != numRefStates_ && referenceStates_.valid())
    {
        FatalError << "different number of reference states and accompanying response functions\n" << abort(FatalError);
    }
    responseFunctions_->readResponseFunctions(dataBaseNames_);
}

fieldNorm& dataBase::fieldN()
{
    if (fieldNorm_.valid())
    {
        return fieldNorm_();
    }
    else
    {
        FatalError << "no field norm set\n" << abort(FatalError);        
    }
}

referenceStates& dataBase::referenceS()
{
    if (referenceStates_.valid())
    {
        return referenceStates_();
    }
    else
    {
        FatalError << "no reference state type set\n" << abort(FatalError); 
    }
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

label dataBase::numRefStates()
{
    return numRefStates_;
}

// fieldListIndex for target field has to be found at solver level via findRefStateListIndex()
/*
label dataBase::findNearestRefState(const volScalarField &field, label fieldListIndex)
{
    scalarList distances(numRefStates_);

    forAll(distances,stateIndex)
    {
        const volScalarField& stateI(referenceS().exportVolScalarField(fieldListIndex,stateIndex));
        distances[stateIndex] = fieldN().fieldsDistance(stateI,field);
    }
    
    scalar minDist = distances[0];
    label minIndex = 0;
    forAll(distances,stateIndex)
    {
        if (distances[stateIndex] < minDist)
        {
            minDist = distances[stateIndex];
            minIndex = stateIndex;
        }
    }

    return minIndex;
}
*/

label dataBase::findNearestRefState(const volVectorField &field, label fieldListIndex, scalar &minDist)
{
    scalarList distances(numRefStates_);

    forAll(distances,stateIndex)
    {
        const volVectorField& stateI(referenceS().exportVolVectorField(fieldListIndex,stateIndex));
        distances[stateIndex] = fieldN().fieldsDistanceConvectiveTerm(stateI,field);
    }
    
    minDist = distances[0];
    label minIndex = 0;
    forAll(distances,stateIndex)
    {
        if (distances[stateIndex] < minDist)
        {
            minDist = distances[stateIndex];
            minIndex = stateIndex;
        }
    }

    return minIndex;
}

labelList& dataBase::faceIDperPatch()
{
    return faceIDperPatch_;
}

labelList& dataBase::patchOwningFace()
{
    return patchOwningFace_;
}

scalar dataBase::predictionTimeStep()
{
    return predictionTimeStep_;
}

void dataBase::setPredictionTimeStep(scalar step)
{
    predictionTimeStep_ = step;
}

}

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Contributing authors:
    Thomas Lichtenegger
    Copyright (C) 2015- Johannes Kepler University, Linz
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
    
    
Class
    Foam::deviationResponse::dataBase

Description
    Handles all data relevant for deviation response CFD calculations.

\*---------------------------------------------------------------------------*/

#ifndef dataBase_H
#define dataBase_H

//#include <vector>

#include "fvCFD.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// forward declarations
class fieldNorm;
class referenceStates;
class responseFunctions;

/*---------------------------------------------------------------------------*\
                           Class dataBase Declaration
\*---------------------------------------------------------------------------*/

class dataBase:
    public regIOobject
{

protected:

    const fvMesh& mesh_;

    IOdictionary dataBaseProperties_;

    fileNameList dataBaseNames_;

    label numDataBases_;

    label numRefStates_;

    scalar predictionTimeStep_;

    autoPtr<fieldNorm> fieldNorm_;

    autoPtr<referenceStates> referenceStates_;

    autoPtr<responseFunctions> responseFunctions_;

    autoPtr<globalIndex> globalCellNumbering_;

    autoPtr<globalIndex> globalBoundaryFaceNumbering_;

    labelList faceIDperPatch_;

    labelList patchOwningFace_;


public:

      // Constructors

    //- Construct from mesh and a list of particles
    dataBase
    (
            const fvMesh& mesh, word dbName = "dataBase"
    );


    //- Destructor
    virtual ~dataBase();

    // public Member Functions

    const fvMesh& mesh() const;

    void init();

    fieldNorm& fieldN();

    referenceStates& referenceS();

    responseFunctions& responseF();

    label localFromGlobalCellID(label cellI);

    label localFromGlobalBoundaryFaceID(label bFaceI);

    label numRefStates();
    
    label findNearestRefState(const volScalarField &field, label fieldListIndex);
    
    label findNearestRefState(const volVectorField &field, label fieldListIndex);

    labelList& faceIDperPatch();

    labelList& patchOwningFace();

    scalar predictionTimeStep();

    void setPredictionTimeStep(scalar);

   // Dummy function for regIOobject
    bool writeData(Ostream&) const
    {
        return true;
    }

};

}

#endif

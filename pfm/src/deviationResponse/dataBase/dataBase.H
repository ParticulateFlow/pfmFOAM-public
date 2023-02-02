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

    wordList dataBaseNames_;

    label numDataBases_;

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
            const fvMesh& mesh
    );


    //- Destructor
    virtual ~dataBase();

    // public Member Functions

    const fvMesh& mesh() const;

    void init();

    referenceStates& referenceS();

    responseFunctions& responseF();

   // Dummy function for regIOobject
    bool writeData(Ostream&) const
    {
        return true;
    }

};

}

#endif

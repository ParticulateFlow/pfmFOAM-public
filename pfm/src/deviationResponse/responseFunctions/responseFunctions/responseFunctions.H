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
    Foam::deviationResponse::responseFunctions::responseFunctions

Description
    Virtual parent class to store response functions.

\*---------------------------------------------------------------------------*/

#ifndef responseFunctions_H
#define responseFunctions_H

#include "dataBase.H"
#include "fvCFD.H"
#include "HashTable.H"
#include "tensorList.H"


namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class responseFunctions Declaration
\*---------------------------------------------------------------------------*/

class responseFunctions
{
protected:

    // Protected data
    dataBase& dataBase_;

    const dictionary& dataBaseProperties_;

    bool verbose_;

    bool readIntegratedResponseFunctions_;

    List<List<labelList>> senderBoundaryFaceIDs_;

    List<List<labelList>> senderCellIDs_;

public:

    friend class recBase;

    //- Runtime type information
    TypeName("responseFunctions");

    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            responseFunctions,
            dictionary,
            (
                const dictionary& dict,
                dataBase& base
            ),
            (dict,base)
        );


    // Constructors

        //- Construct from components
        responseFunctions
        (
            const dictionary& dict,
            dataBase& base
        );


    // Destructor

        virtual ~responseFunctions();


    // Selector

        static autoPtr<responseFunctions> New
        (
            const dictionary& dict,
            dataBase& base
        );


    // Member Functions

    void readSenderIDs(wordList dataBases);

    virtual void readResponseFunctions(wordList dataBases) = 0;

    void readIntegratedResponseFunctions() {readIntegratedResponseFunctions_ = true;}

    labelList &senderCellIDs(label refState, label receiverID); // returns *global* cell IDs sending to cell with *local* receiver ID

    labelList &senderBoundaryFaceIDs(label refState, label receiverID); // returns *global* boundary face IDs sending to cell with *local* receiver ID

    virtual tensorList &Xuu_boundary(label refState, label receiverID) = 0;

    virtual tensorList &Xuu_internal(label refState, label receiverID) = 0;

    virtual tensor Xuu_integrated(label refState, label cellI) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

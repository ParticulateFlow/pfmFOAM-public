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
    Foam::deviationPropagation::deviationPropagators

Description
    Virtual parent class to store deviation propagators.

\*---------------------------------------------------------------------------*/

#ifndef deviationPropagators_H
#define deviationPropagators_H

#include "dataBase.H"
#include "fvCFD.H"
#include "HashTable.H"
#include "tensorList.H"


namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class deviationPropagators Declaration
\*---------------------------------------------------------------------------*/

class deviationPropagators
{
protected:

    // Protected data
    dataBase& dataBase_;

    const dictionary& dataBaseProperties_;

    bool verbose_;

    bool readIntegratedDeviationPropagators_;

    bool deviationPropagatorFormatASCII_;

    List<List<labelList>> senderBoundaryFaceIDs_;

    List<List<labelList>> senderCellIDs_;

public:

    friend class recBase;

    //- Runtime type information
    TypeName("deviationPropagators");

    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            deviationPropagators,
            dictionary,
            (
                const dictionary& dict,
                dataBase& base
            ),
            (dict,base)
        );


    // Constructors

        //- Construct from components
        deviationPropagators
        (
            const dictionary& dict,
            dataBase& base
        );


    // Destructor

        virtual ~deviationPropagators();


    // Selector

        static autoPtr<deviationPropagators> New
        (
            const dictionary& dict,
            dataBase& base
        );


    // Member Functions

    label readSenderIDs(fileNameList dataBases);

    virtual void readDeviationPropagators(fileNameList dataBases) = 0;

    void readIntegratedDeviationPropagators() {readIntegratedDeviationPropagators_ = true;}

    labelList &senderCellIDs(label refState, label receiverID); // returns *global* cell IDs sending to cell with *local* receiver ID

    labelList &senderBoundaryFaceIDs(label refState, label receiverID); // returns *global* boundary face IDs sending to cell with *local* receiver ID

    virtual tensorList &Kuu_boundary(label refState, label receiverID) = 0;

    virtual tensorList &Kuu_internal(label refState, label receiverID) = 0;

    virtual tensor Kuu_integrated(label refState, label cellI) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

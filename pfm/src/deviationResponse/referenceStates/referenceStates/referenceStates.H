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
    Foam::deviationResponse::referenceStates::referenceStates

Description
    Virtual parent class to store reference states.

\*---------------------------------------------------------------------------*/

#ifndef referenceStates_H
#define referenceStates_H

#include "dataBase.H"
#include "fvCFD.H"
#include "HashTable.H"


namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class referenceStates Declaration
\*---------------------------------------------------------------------------*/

class referenceStates
{
protected:

    // Protected data
    dataBase& dataBase_;

    const dictionary& dataBaseProperties_;


    

public:

    friend class recBase;

    //- Runtime type information
    TypeName("referenceStates");

    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            referenceStates,
            dictionary,
            (
                const dictionary& dict,
                dataBase& base
            ),
            (dict,base)
        );


    // Constructors

        //- Construct from components
        referenceStates
        (
            const dictionary& dict,
            dataBase& base
        );


    // Destructor

        virtual ~referenceStates();


    // Selector

        static autoPtr<referenceStates> New
        (
            const dictionary& dict,
            dataBase& base
        );


    // Member Functions

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

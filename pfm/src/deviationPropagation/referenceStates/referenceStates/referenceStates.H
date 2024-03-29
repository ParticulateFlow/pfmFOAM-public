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
    Foam::deviationPropagation::referenceStates::referenceStates

Description
    Virtual parent class to store reference states.

\*---------------------------------------------------------------------------*/

#ifndef referenceStates_H
#define referenceStates_H

#include "dataBase.H"
#include "fvCFD.H"
#include "HashTable.H"
#include "IFstream.H"


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

    scalar stepTolerance_ = 1e-6;

    dictionary propsDict_;

    const wordList volScalarRefStateNames_;

    const wordList volVectorRefStateNames_;


public:

    friend class dataBase;

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
            dataBase& base,
            word type
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

    virtual label readReferenceStates(fileNameList dataBases) = 0;

    label findRefStateListIndex(word fieldType, word fieldName);

    virtual const volScalarField& exportVolScalarField(label, label)  = 0;

    virtual const volVectorField& exportVolVectorField(label, label)  = 0;

    virtual const volScalarField& exportVolScalarEvolvedField(label, label)  = 0;

    virtual const volVectorField& exportVolVectorEvolvedField(label, label)  = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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

#include "error.H"
#include "referenceStates.H"
#include <unistd.h>


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(referenceStates, 0);

defineRunTimeSelectionTable(referenceStates, dictionary);


// * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
referenceStates::referenceStates
(
    const dictionary& dict,
    dataBase& base,
    word type
)
:
    dataBase_(base),
    propsDict_(dict.subDict(type + "Props")),
    volScalarRefStateNames_(propsDict_.lookupOrDefault<wordList>("volScalarRefStates",wordList(0))),
    volVectorRefStateNames_(propsDict_.lookupOrDefault<wordList>("volVectorRefStates",wordList(0)))
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

referenceStates::~referenceStates()
{}

// * * * * * * * * * * * * * public Member Functions  * * * * * * * * * * * * //

label referenceStates::findRefStateListIndex(word fieldType, word fieldName)
{
    if (fieldType == "volScalarField")
    {
        forAll(volScalarRefStateNames_,i)
        {
            if (volScalarRefStateNames_[i] == fieldName)
            {
                return i;
            }
        }
    }
    else if (fieldType == "volVectorField")
    {
        forAll(volVectorRefStateNames_,i)
        {
            if (volVectorRefStateNames_[i] == fieldName)
            {
                return i;
            }
        }
    }
    else
    {
        FatalError << "unknown field type\n" << abort(FatalError);
    }

    FatalError << "could not find field with name " << fieldName << " \n" << abort(FatalError);
    return -1;
}

// * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

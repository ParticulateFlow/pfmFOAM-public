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
#include "responseFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

autoPtr<responseFunctions> responseFunctions::New
(
    const dictionary& dict,
    dataBase& base
)
{
    word responseFunctionsType
    (
        dict.lookup("responseFunctions")
    );

    Info << "Selecting responseFunctions "
         << responseFunctionsType << endl;


    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(responseFunctionsType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "responseFunctions::New(const dictionary&, const spray&) : "
            << endl
            << "    unknown responseFunctionsType type "
            << responseFunctionsType
            << ", constructor not in hash table" << endl << endl
            << "    Valid responseFunctions types are :"
            << endl;
        Info << dictionaryConstructorTablePtr_->toc()
            << abort(FatalError);
    }

    return autoPtr<responseFunctions>(cstrIter()(dict,base));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

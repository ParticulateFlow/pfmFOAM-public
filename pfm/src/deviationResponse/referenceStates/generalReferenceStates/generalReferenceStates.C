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
#include "generalReferenceStates.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(generalReferenceStates, 0);

addToRunTimeSelectionTable
(
    referenceStates,
    generalReferenceStates,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
generalReferenceStates::generalReferenceStates
(
    const dictionary& dict,
    dataBase& base
)
:
    referenceStates(dict,base),
    propsDict_(dict.subDict(typeName + "Props")),
    maxNumRefStates_(propsDict_.lookupOrDefault<label>("maxNumRefStates", 1000)),
    volScalarRefStateList_(volScalarRefStateNames_.size()),
    volVectorRefStateList_(volVectorRefStateNames_.size())
{
    forAll(volScalarRefStateNames_,i)
    {
        volScalarRefStateList_[i].setSize(maxNumRefStates_);
    }

    forAll(volVectorRefStateNames_,i)
    {
        volVectorRefStateList_[i].setSize(maxNumRefStates_);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

generalReferenceStates::~generalReferenceStates()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

label generalReferenceStates::readReferenceStates(wordList dataBases)
{
    int refStates = 0;
    forAll(dataBases,i)
    {
        word dbName = dataBases[i];
        Info << "\nReading reference states of database " << dbName << endl;

        Foam::Time dbTime(fileName(dbName), "", "../system", "../constant", false);
        instantList timeDirs(dbTime.times());
        if (timeDirs.size() == 0)
        {
            FatalError << "database " << dbName << " does not exist or is empty\n" << abort(FatalError);
        }

        for (instantList::iterator it=timeDirs.begin(); it != timeDirs.end(); ++it)
        {
            // set time
            dbTime.setTime(*it, it->value());

            // skip constant
            if (dbTime.timeName() == "constant")
            {
                continue;
            }

            if (verbose_)
            {
                Info << "Reading at t = " << dbTime.timeName() << endl;
            }

            if (refStates >= maxNumRefStates_)
            {
                FatalError << "database contents exceeds maximum number of reference states\n" << abort(FatalError);                
            }

            forAll(volScalarRefStateNames_,j)
            {
                volScalarRefStateList_[j].set
                (
                    refStates,
                    new volScalarField
                    (
                        IOobject
                        (
                            volScalarRefStateNames_[j],
                            dbTime.timePath(),
                            dataBase_.mesh(),
                            IOobject::MUST_READ,
                            IOobject::NO_WRITE
                        ),
                        dataBase_.mesh()
                    )
                );
            }

            forAll(volVectorRefStateNames_,j)
            {
                volVectorRefStateList_[j].set
                (
                    refStates,
                    new volVectorField
                    (
                        IOobject
                        (
                            volVectorRefStateNames_[j],
                            dbTime.timePath(),
                            dataBase_.mesh(),
                            IOobject::MUST_READ,
                            IOobject::NO_WRITE
                        ),
                        dataBase_.mesh()
                    )
                );
            }

            refStates++;
        }
        Info << "Reading reference states of database " << dbName << " done\n" << endl;
    }
    
    forAll(volScalarRefStateNames_,i)
    {
        volScalarRefStateList_[i].resize(refStates);
    }
    
    forAll(volVectorRefStateNames_,i)
    {
        volVectorRefStateList_[i].resize(refStates);
    }

    return refStates;
}

const volScalarField& generalReferenceStates::exportVolScalarField(label fieldListIndex, label refState)
{
    return volScalarRefStateList_[fieldListIndex][refState];
}

const volVectorField& generalReferenceStates::exportVolVectorField(label fieldListIndex, label refState)
{
    return volVectorRefStateList_[fieldListIndex][refState];
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

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
#include "IFstream.H"
#include "responseFunctions.H"
#include <unistd.h>


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(responseFunctions, 0);

defineRunTimeSelectionTable(responseFunctions, dictionary);


// * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
responseFunctions::responseFunctions
(
    const dictionary& dict,
    dataBase& base
)
:
    dataBase_(base),
    dataBaseProperties_(dict),
    verbose_(dict.lookupOrDefault<Switch>("verbose", false)),
    readIntegratedResponseFunctions_(false),
    senderBoundaryFaceIDs_(100000),
    senderCellIDs_(100000)
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

responseFunctions::~responseFunctions()
{}

// * * * * * * * * * * * * * public Member Functions  * * * * * * * * * * * * //

label responseFunctions::readSenderIDs(wordList dataBases)
{
    int refStates = 0;
    forAll(dataBases,i)
    {
        word dbName = dataBases[i];
        Info << "\nReading sender IDs of database " << dbName << endl;

        Foam::Time dbTime(fileName(dbName), "", "../system", "../constant", false);
        instantList timeDirs(dbTime.times());
        if (timeDirs.size() == 0)
        {
            FatalError <<"database " << dbName << " does not exist or is empty\n" << abort(FatalError);
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
            // text files of the form are expected
            // 2
            // (
            // 2(9 1)
            // 3(2 7 4)
            // )
            List<labelList> senderBoundaryFaces;
            IFstream ISBoundaryFaces(dbName+"/"+dbTime.timeName()+"/senderFaces");
            ISBoundaryFaces >> senderBoundaryFaces;
            senderBoundaryFaceIDs_[refStates] = senderBoundaryFaces;

            List<labelList> senderCells;
            IFstream ISCells(dbName+"/"+dbTime.timeName()+"/senderCells");
            ISCells >> senderCells;
            senderCellIDs_[refStates] = senderCells;
            refStates++;
        }
        Info << "Reading sender IDs of database " << dbName <<" done\n" << endl;
    }
    senderBoundaryFaceIDs_.resize(refStates);
    senderCellIDs_.resize(refStates);

    return refStates;
}


labelList& responseFunctions::senderCellIDs(label refState, label receiverID)
{
    // this assumes that the lists of sender IDs are ordered by receiver cell indices; this implies a proper decomposition of these lists for parallel cases
    return senderCellIDs_[refState][receiverID];
}

labelList& responseFunctions::senderBoundaryFaceIDs(label refState, label receiverID)
{
    // this assumes that the lists of sender IDs are ordered by receiver cell indices; this implies a proper decomposition of these lists for parallel cases
    return senderBoundaryFaceIDs_[refState][receiverID];
}
// * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

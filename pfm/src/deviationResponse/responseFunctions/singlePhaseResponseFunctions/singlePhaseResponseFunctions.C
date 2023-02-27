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
#include "singlePhaseResponseFunctions.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(singlePhaseResponseFunctions, 0);

addToRunTimeSelectionTable
(
    responseFunctions,
    singlePhaseResponseFunctions,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
singlePhaseResponseFunctions::singlePhaseResponseFunctions
(
    const dictionary& dict,
    dataBase& base
)
:
    responseFunctions(dict,base),
    propsDict_(dict.subDict(typeName + "Props")),
    Xuu_boundary_(100000),
    Xuu_internal_(100000),
    Xuu_integrated_(100000)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

singlePhaseResponseFunctions::~singlePhaseResponseFunctions()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

label singlePhaseResponseFunctions::readResponseFunctions(wordList dataBases)
{
    int refStates = 0;
    forAll(dataBases,i)
    {
        word dbName = dataBases[i];
        Info << "\nReading response functions of database " << dbName << endl;

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
            // 2((1 2 3 4 5 6 7 8 9) (-1 -2 -3 -4 -5 -6 -7 -8 -9))
            // 3((11 12 13 14 15 16 17 18 19) (21 22 23 24 25 26 27 28 29) (31 32 33 34 35 36 37 38 39))
            // )
            List<tensorList> X_internal;
            IFstream IS_internal(dbName+"/"+dbTime.timeName()+"/X_uu_internal");
            IS_internal >> X_internal;
            Xuu_internal_[refStates] = X_internal;

            List<tensorList> X_boundary;
            IFstream IS_boundary(dbName+"/"+dbTime.timeName()+"/X_uu_boundary");
            IS_boundary >> X_boundary;
            Xuu_boundary_[refStates] = X_boundary;

            if (readIntegratedResponseFunctions_)
            {
                tensorList responsefunction;
                IFstream IS(dbName+"/"+dbTime.timeName()+"/X_uu_integrated");
                IS >> responsefunction;
                Xuu_integrated_[refStates] = responsefunction;
            }

            refStates++;
        }
        Info << "Reading response functions of database " << dbName <<" done\n" << endl;
    }
    Xuu_internal_.resize(refStates);
    Xuu_boundary_.resize(refStates);
    if (readIntegratedResponseFunctions_)
    {
        Xuu_integrated_.resize(refStates);
    }
    else
    {
        Xuu_integrated_.resize(0);
    }

    return refStates;
}

tensorList& singlePhaseResponseFunctions::Xuu_boundary(label refState, label receiverID)
{
    return Xuu_boundary_[refState][receiverID];
}

tensorList& singlePhaseResponseFunctions::Xuu_internal(label refState, label receiverID)
{
    return Xuu_internal_[refState][receiverID];
}

tensor singlePhaseResponseFunctions::Xuu_integrated(label refState, label cellI)
{
    return Xuu_integrated_[refState][cellI];
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

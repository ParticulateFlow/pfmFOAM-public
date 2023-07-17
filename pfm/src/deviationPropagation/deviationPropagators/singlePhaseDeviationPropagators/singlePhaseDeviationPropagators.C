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
#include "singlePhaseDeviationPropagators.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(singlePhaseDeviationPropagators, 0);

addToRunTimeSelectionTable
(
    deviationPropagators,
    singlePhaseDeviationPropagators,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
singlePhaseDeviationPropagators::singlePhaseDeviationPropagators
(
    const dictionary& dict,
    dataBase& base
)
:
    deviationPropagators(dict,base),
    propsDict_(dict.subDict(typeName + "Props")),
    Kuu_boundary_(100000),
    Kuu_internal_(100000),
    Kuu_integrated_(100000)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

singlePhaseDeviationPropagators::~singlePhaseDeviationPropagators()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

label singlePhaseDeviationPropagators::readDeviationPropagators(fileNameList dataBases)
{
    int refStates = 0;
    forAll(dataBases,i)
    {
        fileName dbName = dataBases[i];
        Info << "\nReading deviation propagators of database " << dbName << endl;

        Foam::Time dbTime(dbName, "", "../system", "../constant", false);
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
            List<tensorList> K_internal;
            if (deviationPropagatorFormatASCII_)
            {
                IFstream IS_internal(dbName+"/"+dbTime.timeName()+"/K_uu_internal.txt");
                IS_internal >> K_internal;
            }
            else
            {
                IFstream IS_internal(dbName+"/"+dbTime.timeName()+"/K_uu_internal", IOstream::BINARY);
                IS_internal >> K_internal;
            }
            Kuu_internal_[refStates] = K_internal;

            List<tensorList> K_boundary;
            if (deviationPropagatorFormatASCII_)
            {
                IFstream IS_boundary(dbName+"/"+dbTime.timeName()+"/K_uu_boundary.txt");
                IS_boundary >> K_boundary;
            }
            else
            {
                IFstream IS_boundary(dbName+"/"+dbTime.timeName()+"/K_uu_boundary", IOstream::BINARY);
                IS_boundary >> K_boundary;
            }
            Kuu_boundary_[refStates] = K_boundary;

            if (readIntegratedDeviationPropagators_)
            {
                tensorList K_integrated;
                if (deviationPropagatorFormatASCII_)
                {
                    IFstream IS_integrated(dbName+"/"+dbTime.timeName()+"/K_uu_integrated.txt");
                    IS_integrated >> K_integrated;
                }
                else
                {
                    IFstream IS_integrated(dbName+"/"+dbTime.timeName()+"/K_uu_integrated", IOstream::BINARY);
                    IS_integrated >> K_integrated;
                }
                Kuu_integrated_[refStates] = K_integrated;
            }

            refStates++;
        }
        Info << "Reading deviation propagators of database " << dbName <<" done\n" << endl;
    }
    Kuu_internal_.resize(refStates);
    Kuu_boundary_.resize(refStates);
    if (readIntegratedDeviationPropagators_)
    {
        Kuu_integrated_.resize(refStates);
    }
    else
    {
        Kuu_integrated_.resize(0);
    }

    return refStates;
}

tensorList& singlePhaseDeviationPropagators::Kuu_boundary(label refState, label receiverID)
{
    return Kuu_boundary_[refState][receiverID];
}

tensorList& singlePhaseDeviationPropagators::Kuu_internal(label refState, label receiverID)
{
    return Kuu_internal_[refState][receiverID];
}

tensor singlePhaseDeviationPropagators::Kuu_integrated(label refState, label cellI)
{
    return Kuu_integrated_[refState][cellI];
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

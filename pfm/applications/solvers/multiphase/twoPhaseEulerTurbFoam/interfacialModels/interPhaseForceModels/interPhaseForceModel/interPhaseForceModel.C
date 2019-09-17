/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2018 OpenFOAM Foundation
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

#include "interPhaseForceModel.H"
#include "phasePair.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(interPhaseForceModel, 0);
    defineRunTimeSelectionTable(interPhaseForceModel, dictionary);
}

const Foam::dimensionSet Foam::interPhaseForceModel::dimL(dimless);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interPhaseForceModel::interPhaseForceModel
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject
)
:
    regIOobject
    (
        IOobject
        (
            IOobject::groupName(typeName, pair.name()),
            pair.phase1().mesh().time().timeName(),
            pair.phase1().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            registerObject
        )
    ),
    pair_(pair)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interPhaseForceModel::~interPhaseForceModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::interPhaseForceModel::CpIPh() const
{
    return pair_.dispersed()*Cp();
}

bool Foam::interPhaseForceModel::writeData(Ostream& os) const
{
    return os.good();
}


// ************************************************************************* //

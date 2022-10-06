/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenFOAM Foundation
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

#include "Tenneti.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(Tenneti, 0);
    addToRunTimeSelectionTable(dragModel, Tenneti, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::Tenneti::Tenneti
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject
)
:
    dragModel(dict, pair, registerObject),
    residualRe_("residualRe", dimless, dict.lookup("residualRe"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dragModels::Tenneti::~Tenneti()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::dragModels::Tenneti::CdRe() const
{
    volScalarField alpha1
    (
        max(pair_.dispersed(), pair_.dispersed().residualAlpha())
    );

    volScalarField alpha2
    (
        max(1.0 - pair_.dispersed(), pair_.continuous().residualAlpha())
    );

    volScalarField Res(alpha2*pair_.Re());

    volScalarField F0
    (
        "F0",
         neg(Res - 1000) *(1.0 + 0.15*pow(Res, 0.687))/sqr(alpha2)
       + pos0(Res - 1000)*(0.44/24.0)*Res/sqr(alpha2)
   );

    volScalarField F1
    (
        "F1",
        5.81*(alpha1/sqr(alpha2)) + 0.48*(cbrt(alpha1)/pow3(alpha2))
    );

    volScalarField F2
    (
        "F2",
        Res*alpha2*pow3(alpha1)*(0.95 + (0.61*pow3(alpha1)/sqr(alpha2)))
    );

    return 24.0*alpha2*(F0 + F1 + F2);
}


// ************************************************************************* //

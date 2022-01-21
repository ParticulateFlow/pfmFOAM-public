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

#include "Tang.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(Tang, 0);
    addToRunTimeSelectionTable(dragModel, Tang, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::Tang::Tang
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

Foam::dragModels::Tang::~Tang()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::dragModels::Tang::CdRe() const
{
    volScalarField alpha1
    (
        max(pair_.dispersed(), pair_.dispersed().residualAlpha())
    );

    volScalarField alpha2
    (
        max(scalar(1) - pair_.dispersed(), pair_.continuous().residualAlpha())
    );

    volScalarField Res(alpha2*pair_.Re());

    volScalarField ResLim
    (
        "ReLim",
        max(Res, residualRe_)
    );

    volScalarField F0
    (
        "F0",
        10.0*alpha1/sqr(alpha2) + sqr(alpha2)*(scalar(1.0) + 1.5*sqrt(alpha1))
    );

    volScalarField F1
    (
        "F1",
        Res*(0.11*alpha1*(1.0 + alpha1) - 0.00456/pow(alpha2, 4.0)
      + pow(ResLim, -0.343)*(0.169*alpha2 + 0.0644/pow(alpha2, 4.0)))
    );

    return 24.0*alpha2*(F0 + F1);
}


// ************************************************************************* //

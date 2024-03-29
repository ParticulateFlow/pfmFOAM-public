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

#include "noDriftTemperature.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace driftTemperatureModels
{
    defineTypeNameAndDebug(noDriftTemperature, 0);
    addToRunTimeSelectionTable
    (
        driftTemperatureModel,
        noDriftTemperature,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::driftTemperatureModels::noDriftTemperature::noDriftTemperature
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject
)
:
    driftTemperatureModel(dict, pair, registerObject)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::driftTemperatureModels::noDriftTemperature::
~noDriftTemperature()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::driftTemperatureModels::noDriftTemperature::Tdrift() const
{
    const fvMesh& mesh(this->pair_.phase1().mesh());

    return
        tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "zero",
                    mesh.time().timeName(),
                    mesh
                ),
                mesh,
                dimensionedScalar("zero", dimT, 0)
            )
        );
}


Foam::tmp<Foam::volScalarField>
Foam::driftTemperatureModels::noDriftTemperature::KhTdrift() 
{
    const fvMesh& mesh(this->pair_.phase1().mesh());

    return
        tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "zero",
                    mesh.time().timeName(),
                    mesh
                ),
                mesh,
                dimensionedScalar("zero", dimG, 0)
            )
        );
}




// ************************************************************************* //

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

#include "noneRegularization.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace ADMdispersedModels
{
    defineTypeNameAndDebug(noneRegularization, 0);
    addToRunTimeSelectionTable(regularizationModel, noneRegularization, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ADMdispersedModels::noneRegularization::noneRegularization
(
    const dictionary& dict,
    const volScalarField& alpha
)
:
    regularizationModel(dict,alpha)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ADMdispersedModels::noneRegularization::~noneRegularization()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField> Foam::ADMdispersedModels::noneRegularization::regTerm
(
     const volScalarField& alpha,
     const volScalarField& rho,
     const volScalarField& k,
     const volVectorField& U,
     const volVectorField& Ustar
) const
{
    return dimensionedVector
    (
        "0",
        dimensionSet(1, -2, -2, 0, 0, 0, 0),
        vector(0.0,0.0,0.0)
    )*alpha;
}


// ************************************************************************* //

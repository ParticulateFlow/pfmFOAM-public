/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "ParmentierADMFilter.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ParmentierADMFilter, 0);
    addToRunTimeSelectionTable(LESfilter, ParmentierADMFilter, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ParmentierADMFilter::ParmentierADMFilter
(
    const fvMesh& mesh
)
:
    LESfilter(mesh)
{}


Foam::ParmentierADMFilter::ParmentierADMFilter(const fvMesh& mesh, const dictionary&)
:
    LESfilter(mesh)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ParmentierADMFilter::read(const dictionary&)
{}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::ParmentierADMFilter::operator()
(
    const tmp<volScalarField>& unFilteredField
) const
{
    correctBoundaryConditions(unFilteredField);

    tmp<volScalarField> filteredField =
          0.5*unFilteredField()
        + 0.5*fvc::surfaceSum
          (
              fvc::interpolate(unFilteredField())
          )
        / fvc::surfaceSum(mesh().magSf()/mesh().magSf());

    unFilteredField.clear();

    return filteredField;
}


Foam::tmp<Foam::volVectorField> Foam::ParmentierADMFilter::operator()
(
    const tmp<volVectorField>& unFilteredField
) const
{
    correctBoundaryConditions(unFilteredField);
    
    tmp<volVectorField> filteredField =
          0.5*unFilteredField()
        + 0.5*fvc::surfaceSum
          (
              fvc::interpolate(unFilteredField())
          )
        / fvc::surfaceSum(mesh().magSf()/mesh().magSf());


    return filteredField;
}


Foam::tmp<Foam::volSymmTensorField> Foam::ParmentierADMFilter::operator()
(
    const tmp<volSymmTensorField>& unFilteredField
) const
{
    correctBoundaryConditions(unFilteredField);
    
    tmp<volSymmTensorField> filteredField =
          0.5*unFilteredField()
        + 0.5*fvc::surfaceSum
          (
              fvc::interpolate(unFilteredField())
          )
        / fvc::surfaceSum(mesh().magSf()/mesh().magSf());

    unFilteredField.clear();

    return filteredField;
}


Foam::tmp<Foam::volTensorField> Foam::ParmentierADMFilter::operator()
(
    const tmp<volTensorField>& unFilteredField
) const
{
    correctBoundaryConditions(unFilteredField);

    tmp<volTensorField> filteredField =
          0.5*unFilteredField()
        + 0.5*fvc::surfaceSum
          (
              fvc::interpolate(unFilteredField())
          )
        / fvc::surfaceSum(mesh().magSf()/mesh().magSf());

    unFilteredField.clear();

    return filteredField;
}


// ************************************************************************* //

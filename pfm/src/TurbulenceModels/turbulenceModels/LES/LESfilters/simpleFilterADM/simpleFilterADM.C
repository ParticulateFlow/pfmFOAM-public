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

#include "simpleFilterADM.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(simpleFilterADM, 0);
    addToRunTimeSelectionTable(LESfilter, simpleFilterADM, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::simpleFilterADM::simpleFilterADM
(
    const fvMesh& mesh
)
:
    LESfilter(mesh)
{}


Foam::simpleFilterADM::simpleFilterADM(const fvMesh& mesh, const dictionary&)
:
    LESfilter(mesh)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::simpleFilterADM::read(const dictionary&)
{}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::simpleFilterADM::operator()
(
    const tmp<volScalarField>& unFilteredField
) const
{
    correctBoundaryConditions(unFilteredField);

    tmp<volScalarField> filteredField = 0.5*unFilteredField()
    + (1.0/2.0)*fvc::surfaceSum
    (
        mesh().magSf()*fvc::interpolate(unFilteredField())
    )/fvc::surfaceSum(mesh().magSf());

    unFilteredField.clear();

    return filteredField;
}


Foam::tmp<Foam::volVectorField> Foam::simpleFilterADM::operator()
(
    const tmp<volVectorField>& unFilteredField
) const
{
    correctBoundaryConditions(unFilteredField);

    tmp<volVectorField> filteredField = 0.5*unFilteredField()
    + (1.0/2.0)*fvc::surfaceSum
    (
        mesh().magSf()*fvc::interpolate(unFilteredField())
    )/fvc::surfaceSum(mesh().magSf());
    unFilteredField.clear();

    return filteredField;
}


Foam::tmp<Foam::volSymmTensorField> Foam::simpleFilterADM::operator()
(
    const tmp<volSymmTensorField>& unFilteredField
) const
{
    correctBoundaryConditions(unFilteredField);

    tmp<volSymmTensorField> filteredField = 0.5*unFilteredField()
    + (1.0/2.0)*fvc::surfaceSum
    (
        mesh().magSf()*fvc::interpolate(unFilteredField())
    )/fvc::surfaceSum(mesh().magSf());

    unFilteredField.clear();

    return filteredField;
}


Foam::tmp<Foam::volTensorField> Foam::simpleFilterADM::operator()
(
    const tmp<volTensorField>& unFilteredField
) const
{
    correctBoundaryConditions(unFilteredField);

    tmp<volTensorField> filteredField = 0.5*unFilteredField()
    + (1.0/2.0)*fvc::surfaceSum
    (
        mesh().magSf()*fvc::interpolate(unFilteredField())
    )/fvc::surfaceSum(mesh().magSf());

    unFilteredField.clear();

    return filteredField;
}


// ************************************************************************* //

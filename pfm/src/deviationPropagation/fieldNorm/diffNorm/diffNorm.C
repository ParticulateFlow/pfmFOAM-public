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
#include "diffNorm.H"
#include "referenceStates.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(diffNorm, 0);

addToRunTimeSelectionTable
(
    fieldNorm,
    diffNorm,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
diffNorm::diffNorm
(
    const dictionary& dict,
    dataBase& base,
    word type
)
:
    fieldNorm(dict, base),
    propsDict_(dict.subDict(type + "Props")),
    verbose_(propsDict_.lookupOrDefault<bool>("verbose", false)),
    fieldName_(propsDict_.lookup("fieldName")),
    domainVolume_(0.0)
{
    domainVolume_ = gSum(dataBase_.mesh().V());
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

diffNorm::~diffNorm()
{}

// * * * * * * * * * * * * * protected Member Functions  * * * * * * * * * * * * //

scalar diffNorm::fieldsDistance(const volScalarField &field1, const volScalarField &field2, double normalization)
{
    volScalarField diffField(field1-field2);
    scalar integrand = 0.0;
    scalar sum = 0.0;
    forAll(diffField,cellI)
    {
        integrand = sqr(diffField[cellI]);
        if (integrand > VSMALL)
        {
            sum += Foam::sqrt(integrand) * dataBase_.mesh().V()[cellI];
        }
    }
    reduce(sum, sumOp<scalar>());
    // if no normalization value is provided, normalize with domain volume
    if (normalization < -0.5) return sum/domainVolume_;
    else if (fabs(normalization) < VSMALL) FatalError << "attempted normalization with 0 or very small value\n" << abort(FatalError);
    return sum / normalization;
}

scalar diffNorm::fieldsDistance(const volVectorField &field1, const volVectorField &field2, double normalization)
{
    volVectorField diffField(field1-field2);
    scalar integrand = 0.0;
    scalar sum = 0.0;
    forAll(diffField,cellI)
    {
        integrand = magSqr(diffField[cellI]);
        if (integrand > VSMALL)
        {
            sum += Foam::sqrt(integrand) * dataBase_.mesh().V()[cellI];
        }
    }
    reduce(sum, sumOp<scalar>());
    // if no normalization value is provided, normalize with domain volume
    if (normalization < -0.5) return sum/domainVolume_;
    else if (fabs(normalization) < VSMALL) FatalError << "attempted normalization with 0 or very small value\n" << abort(FatalError);
    return sum / normalization;
}

scalar diffNorm::fieldsDistance(const volTensorField &field1, const volTensorField &field2, double normalization)
{
    volTensorField diffField(field1-field2);
    scalar integrand = 0.0;
    scalar sum = 0.0;
    forAll(diffField,cellI)
    {
        integrand = magSqr(diffField[cellI]);
        if (integrand > VSMALL)
        {
            sum += Foam::sqrt(integrand) * dataBase_.mesh().V()[cellI];
        }
    }
    reduce(sum, sumOp<scalar>());
    // if no normalization value is provided, normalize with domain volume
    if (normalization < -0.5) return sum/domainVolume_;
    else if (fabs(normalization) < VSMALL) FatalError << "attempted normalization with 0 or very small value\n" << abort(FatalError);
    return sum / normalization;
}

scalar diffNorm::fieldsDistanceConvectiveTerm(const volVectorField &field1, const volVectorField &field2, double normalization)
{
    volVectorField diffField(field1-field2);
    surfaceScalarField diffFieldSurface(linearInterpolate(diffField) & dataBase_.mesh().Sf());
    volVectorField convectiveTermError(fvc::div(diffFieldSurface,diffField,"convectiveTermError"));
    scalar integrand = 0.0;
    scalar sum = 0.0;
    forAll(diffField,cellI)
    {
        integrand = magSqr(convectiveTermError[cellI]);
        if (integrand > VSMALL)
        {
            sum += Foam::sqrt(integrand) * dataBase_.mesh().V()[cellI];
        }
    }
    reduce(sum, sumOp<scalar>());
    // if no normalization value is provided, normalize with domain volume
    if (normalization < -0.5) return sum/domainVolume_;
    else if (fabs(normalization) < VSMALL) FatalError << "attempted normalization with 0 or very small value\n" << abort(FatalError);
    return sum / normalization;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

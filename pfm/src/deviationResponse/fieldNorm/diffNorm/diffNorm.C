/*---------------------------------------------------------------------------*\
    CFDEMcoupling academic - Open Source CFD-DEM coupling

    Contributing authors:
    Thomas Lichtenegger
    Copyright (C) 2015- Johannes Kepler University, Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling academic.

    CFDEMcoupling academic is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CFDEMcoupling academic is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling academic.  If not, see <http://www.gnu.org/licenses/>.
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
    fieldName_(propsDict_.lookup("fieldName"))
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

diffNorm::~diffNorm()
{}

// * * * * * * * * * * * * * protected Member Functions  * * * * * * * * * * * * //

scalar diffNorm::fieldsDistance(const volScalarField &field1, const volScalarField &field2)
{
    volScalarField diffField(field1-field2);
    return fvc::domainIntegrate(sqr(diffField)).value();
}

scalar diffNorm::fieldsDistance(const volVectorField &field1, const volVectorField &field2)
{
    volVectorField diffField(field1-field2);
    return fvc::domainIntegrate(magSqr(diffField)).value();
}

scalar diffNorm::fieldsDistance(const volTensorField &field1, const volTensorField &field2)
{
    volTensorField diffField(field1-field2);
    return fvc::domainIntegrate(magSqr(diffField)).value();
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

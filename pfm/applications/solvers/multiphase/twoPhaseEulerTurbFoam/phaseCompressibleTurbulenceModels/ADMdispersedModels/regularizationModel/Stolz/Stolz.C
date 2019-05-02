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

#include "Stolz.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"
#include "calculatedFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace ADMdispersedModels
{
namespace regularizationModels
{
    defineTypeNameAndDebug(Stolz, 0);
    addToRunTimeSelectionTable(regularizationModel, Stolz, dictionary);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ADMdispersedModels::regularizationModels::Stolz::Stolz
(
    const dictionary& dict,
    const volScalarField& alpha
)
:
    regularizationModel(dict,alpha),
    coeffDict_(dict.optionalSubDict(typeName + "Coeffs")),
    smagConst_("smagConst", dimless, coeffDict_),
    lengthConst_("lengthConst", dimless, coeffDict_),
    filterPtr_(LESfilter::New(alpha.mesh(), coeffDict_)),
    filter_(filterPtr_())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ADMdispersedModels::regularizationModels::Stolz::~Stolz()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField>
Foam::ADMdispersedModels::regularizationModels::Stolz::regTerm
(
     const volScalarField& alpha,
     const volScalarField& rho,
     const volScalarField& k,
     const volVectorField& U,
     const volVectorField& Ustar
) const
{
    volScalarField V
    (
        IOobject
        (
            "V",
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
         ),
        U.mesh(),
        dimensionedScalar("small", dimLength*dimLength*dimLength, 1.e-2),
        calculatedFvPatchScalarField::typeName
    );
    V.ref() = U.mesh().V();
    volScalarField delta(pow(V,1.0/3.0));
    
    volScalarField tau_1 = sqrt(k)/(lengthConst_*smagConst_*delta);
    
    return alpha*rho*(U - filter_(Ustar))*tau_1;
}

bool Foam::ADMdispersedModels::regularizationModels::
Stolz::read()
{
    coeffDict_ <<= dict_.optionalSubDict(typeName + "Coeffs");
    
    smagConst_.read(coeffDict_);
    lengthConst_.read(coeffDict_);
   
    return true;
}


// ************************************************************************* //

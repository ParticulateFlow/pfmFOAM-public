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

#include "SchneiderbauerEtAlPressure.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace granularPressureModels
{
    defineTypeNameAndDebug(SchneiderbauerEtAl, 0);

    addToRunTimeSelectionTable
    (
        granularPressureModel,
        SchneiderbauerEtAl,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::granularPressureModels::SchneiderbauerEtAl::SchneiderbauerEtAl
(
    const dictionary& dict
)
:
    granularPressureModel(dict),
    coeffDict_(dict.optionalSubDict(typeName + "Coeffs")),
    L_("L", dimensionSet(0, 1, 0, 0, 0), coeffDict_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::granularPressureModels::SchneiderbauerEtAl::~SchneiderbauerEtAl()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::granularPressureModels::SchneiderbauerEtAl::granularPressureCoeff
(
    const volScalarField& alpha1,
    const volScalarField& g0,
    const volScalarField& rho1,
    const volScalarField& da,
    const dimensionedScalar& e
) const
{
    volScalarField lambda
    (
         scalar(1.0) + da/(6.0*sqrt(2.0)*(alpha1 + scalar(1.0e-7))*L_)
    );
    
    return rho1*alpha1*(1.0/lambda + 2.0*(1.0 + e)*alpha1*g0);
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::granularPressureModels::SchneiderbauerEtAl::
granularPressureCoeffPrime
(
    const volScalarField& alpha1,
    const volScalarField& g0,
    const volScalarField& g0prime,
    const volScalarField& rho1,
    const volScalarField& da,
    const dimensionedScalar& e
) const
{
    volScalarField lambda
    (
         scalar(1.0) + da/(6.0*sqrt(2.0)*(alpha1 + scalar(1.0e-7))*L_)
    );
    volScalarField lambdaP
    (
        - da/(6.0*sqrt(2.0)*sqr(alpha1 + scalar(1.0e-7))*L_)
    );
    return rho1
        *(
            2.0*(1.0 + e)*alpha1*g0
          + 1.0/lambda
          + alpha1
           *(
                2.0*(1.0 + e)*g0
              + 2.0*(1.0 + e)*alpha1*g0prime
              - lambdaP/sqr(lambdaP)
            )
         );
}


// ************************************************************************* //

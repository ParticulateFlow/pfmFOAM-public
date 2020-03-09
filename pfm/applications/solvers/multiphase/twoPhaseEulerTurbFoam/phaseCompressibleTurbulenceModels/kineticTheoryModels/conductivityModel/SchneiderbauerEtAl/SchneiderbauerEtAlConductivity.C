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

#include "SchneiderbauerEtAlConductivity.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace conductivityModels
{
    defineTypeNameAndDebug(SchneiderbauerEtAl, 0);

    addToRunTimeSelectionTable
    (
        conductivityModel,
        SchneiderbauerEtAl,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::conductivityModels::SchneiderbauerEtAl::SchneiderbauerEtAl
(
    const dictionary& dict
)
:
    conductivityModel(dict),
    coeffDict_(dict.optionalSubDict(typeName + "Coeffs")),
    L_("L", dimensionSet(0, 1, 0, 0, 0), coeffDict_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::conductivityModels::SchneiderbauerEtAl::
~SchneiderbauerEtAl()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::conductivityModels::SchneiderbauerEtAl::kappa
(
    const volScalarField& alpha1,
    const volScalarField& Theta,
    const volScalarField& g0,
    const volScalarField& rho1,
    const volScalarField& da,
    const dimensionedScalar& e
) const
{
    const scalar Pi = constant::mathematical::pi;
    const scalar sqrtPi = sqrt(Pi);
    
    const fvMesh& mesh = alpha1.mesh();
    
    const volScalarField& Kd = mesh.lookupObject<volScalarField>("Kd");

    volScalarField lambda
    (
        scalar(1) + da/(6.0*sqrt(2.0)*(alpha1 + scalar(1.0e-7)))/L_
     );
    
    dimensionedScalar eta(0.5*(1 + e));
    
    volScalarField kappa
    (
        1.5625*rho1*da*sqrtPi*sqrt(Theta)
       /(eta*(scalar(41.0) - 33.0*eta))
    );
    
    volScalarField kappaStar
    (
        kappa
       /(
            scalar(1.0)
          + 6.0*Kd*kappa
           /(5.0*sqr(alpha1*rho1)*g0*Theta)
         )
    );

    return
     kappaStar/g0
    *(
         (
            1.0/lambda
          + 2.4*eta*alpha1*g0
         )
        *(
            scalar(1.0)
          + 2.4*sqr(eta)*(4.0*eta - scalar(3.0)*alpha1*g0)
         )
       + (2.56/Pi)*(scalar(41.0) - 33.0*eta)*sqr(eta*alpha1*g0)
    );
}


bool Foam::kineticTheoryModels::conductivityModels::SchneiderbauerEtAl::read()
{
    coeffDict_ <<= dict_.optionalSubDict(typeName + "Coeffs");

    L_.readIfPresent(coeffDict_);

    return true;
}


// ************************************************************************* //

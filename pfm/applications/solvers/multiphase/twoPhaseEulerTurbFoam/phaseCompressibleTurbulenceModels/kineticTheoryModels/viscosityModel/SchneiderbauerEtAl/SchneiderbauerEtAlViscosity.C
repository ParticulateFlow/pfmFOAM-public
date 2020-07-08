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

#include "SchneiderbauerEtAlViscosity.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace viscosityModels
{
    defineTypeNameAndDebug(SchneiderbauerEtAl, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        SchneiderbauerEtAl,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::viscosityModels::SchneiderbauerEtAl::SchneiderbauerEtAl
(
    const dictionary& dict
)
:
    viscosityModel(dict),
    coeffDict_(dict.optionalSubDict(typeName + "Coeffs")),
    L_("L", dimensionSet(0, 1, 0, 0, 0), coeffDict_),
    alphaPre_("alpha", dimensionSet(0, 0, 0, 0, 0), coeffDict_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::viscosityModels::SchneiderbauerEtAl::~SchneiderbauerEtAl()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::viscosityModels::SchneiderbauerEtAl::nu
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

    volScalarField mu((5.0/96.0)*rho1*da*sqrtPi*sqrt(Theta));
    volScalarField muStar
    (
        mu
      /(
           scalar(1.0)
         + 2.0*Kd*mu
          /(sqr((alpha1+scalar(1.0e-7))*rho1)*g0*Theta)
       )
    );
    volScalarField muB(51.2*mu*sqr(alpha1)*g0/Pi);

    volScalarField lambda
    (
        scalar(1) + da/(6.0*sqrt(2.0)*(alpha1 + scalar(1.0e-7)))/L_
    );
    
    dimensionedScalar eta(0.5*(1 + e));
    
    return
     (scalar(2.0) + alphaPre_)/(3.0*rho1)
    *(
        muStar
      /(g0*eta*(scalar(2.0) - eta))
      *(
            1.0/lambda
          + 1.6*alpha1*eta*g0
       )
      *(
           scalar(1.0)
         + 1.6*eta*(3.0*eta - scalar(2.0))*alpha1*g0
       )
      + 0.6*eta*muB
    );
}


bool Foam::kineticTheoryModels::viscosityModels::SchneiderbauerEtAl::read()
{
    coeffDict_ <<= dict_.optionalSubDict(typeName + "Coeffs");

    L_.readIfPresent(coeffDict_);
    alphaPre_.readIfPresent(coeffDict_);

    return true;
}


// ************************************************************************* //

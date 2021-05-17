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

Class
    Foam::kineticTheoryModels::viscosityModels::SchneiderbauerEtAl

Description

SourceFiles
    SchneiderbauerEtAlViscosity.C

\*---------------------------------------------------------------------------*/

#ifndef SchneiderbauerEtAlViscosity_H
#define SchneiderbauerEtAlViscosity_H

#include "viscosityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace viscosityModels
{

/*---------------------------------------------------------------------------*\
                           Class SchneiderbauerEtAl Declaration
\*---------------------------------------------------------------------------*/

class SchneiderbauerEtAl
:
    public viscosityModel
{
    // Private data

        dictionary coeffDict_;

        //- Characteristic length of geometry
        dimensionedScalar L_;
    
        //- prefactor
        dimensionedScalar alphaPre_;


public:

    //- Runtime type information
    TypeName("SchneiderbauerEtAl");


    // Constructors

        //- Construct from components
        SchneiderbauerEtAl(const dictionary& dict);


    //- Destructor
    virtual ~SchneiderbauerEtAl();


    // Member functions

        tmp<volScalarField> nu
        (
            const volScalarField& alpha1,
            const volScalarField& Theta,
            const volScalarField& g0,
            const volScalarField& rho1,
            const volScalarField& da,
            const dimensionedScalar& e
        ) const;

        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace viscosityModels
} // End namespace kineticTheoryModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

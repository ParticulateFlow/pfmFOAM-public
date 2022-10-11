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
    Foam::turbGenPope

Description
    Generate a turbulent velocity field conforming to a given
    energy spectrum and being divergence free.

SourceFiles
    turbGenPope.C
    turbGenPopeIO.C

\*---------------------------------------------------------------------------*/

#ifndef turbGenPope_H
#define turbGenPope_H

#include "Random.H"
#include "primitiveFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Kmesh;

/*---------------------------------------------------------------------------*\
                           Class turbGenPope Declaration
\*---------------------------------------------------------------------------*/

class turbGenPope
{
    // Private data

        const Kmesh& K;
        const scalar Cp;
        const scalar p0;
        const scalar beta;
        const scalar nu;
        const scalar cL;
        const scalar cEta;
        const scalar e;
        const scalar epsilon;
        Random RanGen;


public:

    // Constructors

        //- Construct from components
        turbGenPope
        (
             const Kmesh& k,
             const scalar CP,
             const scalar P0,
             const scalar BETA,
             const scalar NU,
             const scalar CL,
             const scalar CETA,
             const scalar E,
             const scalar EPSILON
        );


    // Member Functions

        //- Generate and return a velocity field
        vectorField U();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

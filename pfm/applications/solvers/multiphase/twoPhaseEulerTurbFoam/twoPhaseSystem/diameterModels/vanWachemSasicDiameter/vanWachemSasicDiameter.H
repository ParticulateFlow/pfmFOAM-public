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
    Foam::diameterModels::vanWachemSasic

Description
    Force balance to predict agglomerate size of cohesive particles
    \verbatim .
       Wachem, B.G.M. van, S. Sasic. Derivation, Simulation and Validation of a Cohesive Particle Flow CFD Model. AIChE J., 2008, 54 (1), 9–19. DOI: 10.1002/aic.
    \endverbatim

SourceFiles
    vanWachemSasic.C

\*---------------------------------------------------------------------------*/

#ifndef vanWachemSasicDiameter_H
#define vanWachemSasicDiameter_H

#include "diameterModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{

/*---------------------------------------------------------------------------*\
                           Class vanWachemSasic Declaration
\*---------------------------------------------------------------------------*/

class vanWachemSasic
:
    public diameterModel
{
    // Private data

        //- The vanWachemSasic diameter of the phase
        dimensionedScalar d0_;
        dimensionedScalar H_;
        dimensionedScalar delta_;
        dimensionedScalar rhoB_;
        dimensionedScalar g_;
        dimensionedScalar nu_;
        dimensionedScalar E_;
    
        //- The Sauter-mean diameter of the phase
       volScalarField d_;


public:

    //- Runtime type information
    TypeName("vanWachemSasic");


    // Constructors

        //- Construct from components
        vanWachemSasic
        (
            const dictionary& diameterProperties,
            const phaseModel& phase
        );


    //- Destructor
    virtual ~vanWachemSasic();


    // Member Functions

        //- Return the diameter as a field
        virtual tmp<volScalarField> d() const
        {
            return d_;
        }

        //- Correct the diameter field
        virtual void correct();

        //- Read diameterProperties dictionary
        virtual bool read(const dictionary& diameterProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

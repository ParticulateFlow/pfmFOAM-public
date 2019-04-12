/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2018 OpenFOAM Foundation
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
    Foam::driftVelocityModels::noDriftVelocity

Description

SourceFiles
    noDriftVelocity.C

\*---------------------------------------------------------------------------*/

#ifndef noDriftVelocity_H
#define noDriftVelocity_H

#include "driftVelocityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace driftVelocityModels
{

/*---------------------------------------------------------------------------*\
                      Class noDriftVelocity Declaration
\*---------------------------------------------------------------------------*/

class noDriftVelocity
:
    public driftVelocityModel
{
public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from a dictionary and a phase pair
        noDriftVelocity
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~noDriftVelocity();


    // Member Functions

        //- Turbulent drift velocity
        virtual tmp<volVectorField> udrift() const;
        //- Turbulent dispersion force
        virtual tmp<volVectorField> KdUdrift() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace driftVelocityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

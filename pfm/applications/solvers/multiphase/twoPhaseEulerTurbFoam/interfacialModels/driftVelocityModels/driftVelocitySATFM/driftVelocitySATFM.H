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
    Foam::driftVelocityModels::driftVelocitySATFM
 
(c) Simon Schneiderbauer 2019
    Christian Doppler Laboratory for Multiscale Modeling of Multiphase Processes
    Johannes Kepler University, Linz, Austria

Description
    Turbulent drift velocity model of Schneiderbauer & Saeedipour

    References:
    \verbatim
        Schneiderbauer, S., & Saeedipour, M. (2018). Approximate deconvolution model for the simulation of turbulent gas-solid flows: An a-priori analysis. Physics of Fluids, 30(2), 023301.
    \endverbatim

SourceFiles
    driftVelocitySATFM.C

\*---------------------------------------------------------------------------*/

#ifndef driftVelocitySATFM_H
#define driftVelocitySATFM_H

#include "driftVelocityModel.H"
#include "LESfilter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace driftVelocityModels
{

/*---------------------------------------------------------------------------*\
                             Class driftVelocitySATFM Declaration
\*---------------------------------------------------------------------------*/

class driftVelocitySATFM
:
    public driftVelocityModel
{
    
    //- Residual phase fraction
    const dimensionedScalar residualAlpha_;

public:

    
    //- Runtime type information
    TypeName("driftVelocitySATFM");


    // Constructors

        //- Construct from a dictionary and a phase pair
        driftVelocitySATFM
        (
            const dictionary& dict,
            const phasePair& pair,
            const bool registerObject
        );


    //- Destructor
    virtual ~driftVelocitySATFM();


    // Member Functions

        //- Turbulent drift velocity
        virtual tmp<volVectorField> udrift() const;
    
        //- compute drift velocity
        // virtual void computeDriftVel();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace driftVelocitySATFM
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

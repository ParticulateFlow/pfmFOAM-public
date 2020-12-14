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
    Foam::driftVelocityModels::driftVelocitySATFMreg
 
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
    driftVelocitySATFMreg.C

\*---------------------------------------------------------------------------*/

#ifndef driftVelocitySATFMreg_H
#define driftVelocitySATFMreg_H

#include "driftVelocityModel.H"
#include "LESfilter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace driftVelocityModels
{

/*---------------------------------------------------------------------------*\
                             Class driftVelocitySATFMreg Declaration
\*---------------------------------------------------------------------------*/

class driftVelocitySATFMreg
:
    public driftVelocityModel
{
    
    //- Residual phase fraction
    const dimensionedScalar residualAlpha_;

public:

    
    //- Runtime type information
    TypeName("driftVelocitySATFMreg");


    // Constructors

        //- Construct from a dictionary and a phase pair
        driftVelocitySATFMreg
        (
            const dictionary& dict,
            const phasePair& pair,
            const bool registerObject
        );


    //- Destructor
    virtual ~driftVelocitySATFMreg();


    // Member Functions

        //- Turbulent drift velocity
        virtual tmp<volVectorField> udrift() const;
    
        //- compute drift velocity
        // virtual void computeDriftVel();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace driftVelocitySATFMreg
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

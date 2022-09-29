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
    Foam::driftTemperatureModels::driftTemperatureSATFM
 
(c) 2022
    Christian Doppler Laboratory for Multiscale Modeling of Multiphase Processes
    Johannes Kepler University, Linz, Austria

Description
    Drift Temperature model of Rauchenzauner & Schneiderbauer

    References:
    \verbatim
        Rauchenzauner, S. and Schneiderbauer, S. (2020). A Dynamic Spatially-Averaged Two-Fluid Model for Heat Transport in Moderately Dense Gas-Particle Flows. Physics of Fluids, 32, 063307.
        Rauchenzauner, S. and Schneiderbauer, S. (2022). Validation study of a Spatially-Averaged Two-Fluid Model for heat transport in gas-particle flows. Int. J. Heat Mass Transf., 198, 123382.
    \endverbatim

SourceFiles
    driftTemperatureSATFM.C

\*---------------------------------------------------------------------------*/

#ifndef driftTemperatureSATFM_H
#define driftTemperatureSATFM_H

#include "driftTemperatureModel.H"
#include "LESfilter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace driftTemperatureModels
{

/*---------------------------------------------------------------------------*\
                             Class driftTemperatureSATFM Declaration
\*---------------------------------------------------------------------------*/

class driftTemperatureSATFM
:
    public driftTemperatureModel
{
    
    //- Residual phase fraction
    const dimensionedScalar residualAlpha_;

public:

    
    //- Runtime type information
    TypeName("driftTemperatureSATFM");


    // Constructors

        //- Construct from a dictionary and a phase pair
        driftTemperatureSATFM
        (
            const dictionary& dict,
            const phasePair& pair,
            const bool registerObject
        );


    //- Destructor
    virtual ~driftTemperatureSATFM();


    // Member Functions

        //- Turbulent drift Temperature
        virtual tmp<volScalarField> Tdrift() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace driftTemperatureSATFM
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

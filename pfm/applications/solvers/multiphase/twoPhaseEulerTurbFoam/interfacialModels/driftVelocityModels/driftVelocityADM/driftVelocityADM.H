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
    Foam::driftVelocityModels::driftVelocityADM
 
(c) Simon Schneiderbauer 2019
    Christian Doppler Laboratory for Multiscale Modeling of Multiphase Processes
    Johannes Kepler University, Linz, Austria

Description
    Turbulent drift velocity model of Schneiderbauer & Saeedipour

    References:
    \verbatim
        Schneiderbauer, S., & Saeedipour, M. (2019). "Numerical simulation of turbulent gas-solid flow using an approximate deconvolution model." International Journal of Multiphase Flow, 114, 287–302.
    \endverbatim

SourceFiles
    driftVelocityADM.C

\*---------------------------------------------------------------------------*/

#ifndef driftVelocityADM_H
#define driftVelocityADM_H

#include "driftVelocityModel.H"
#include "LESfilter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace driftVelocityModels
{

/*---------------------------------------------------------------------------*\
                             Class driftVelocityADM Declaration
\*---------------------------------------------------------------------------*/

class driftVelocityADM
:
    public driftVelocityModel
{
    
    //- Residual phase fraction
    const dimensionedScalar residualAlpha_;
    
    //- max phase fraction
    const dimensionedScalar alphaMax_;
    
protected:
    autoPtr<LESfilter> filterPtr_;
    LESfilter& filter_;

public:

    
    //- Runtime type information
    TypeName("driftVelocityADM");


    // Constructors

        //- Construct from a dictionary and a phase pair
        driftVelocityADM
        (
            const dictionary& dict,
            const phasePair& pair,
            const bool registerObject
        );


    //- Destructor
    virtual ~driftVelocityADM();


    // Member Functions

        //- Turbulent drift velocity
        virtual tmp<volVectorField> udrift() const;
    
        //- compute drift velocity
        // virtual void computeDriftVel();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace driftVelocityADM
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

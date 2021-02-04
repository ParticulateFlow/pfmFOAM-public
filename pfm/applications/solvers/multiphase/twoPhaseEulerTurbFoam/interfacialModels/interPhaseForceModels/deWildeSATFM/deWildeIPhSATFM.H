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
    Foam::interPhaseForceModels::deWildeIPhSATFM

(c) Simon Schneiderbauer 2019
    Christian Doppler Laboratory for Multiscale Modeling of Multiphase Processes
    Johannes Kepler University, Linz, Austria
 
Description
    Virtual mass model of de Wilde.

    Reference:
    \verbatim
        De Wilde, J. (2005). Reformulating and quantifying the generalized added mass in filtered gas-solid flow models. Physics of Fluids, 17, 113304.
        Rauchenzauner & Schneiderbauer (2019). IJMF.
    \endverbatim

SourceFiles
    deWildeIPhSATFM.C

\*---------------------------------------------------------------------------*/

#ifndef deWildeIPhSATFM_H
#define deWildeIPhSATFM_H

#include "interPhaseForceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace interPhaseForceModels
{

/*---------------------------------------------------------------------------*\
                           Class deWildeIPhSATFM Declaration
\*---------------------------------------------------------------------------*/

class deWildeIPhSATFM
:
    public interPhaseForceModel
{
    //- Residual phase fraction
    const dimensionedScalar residualAlpha_;

public:

    //- Runtime type information
    TypeName("deWildeIPhSATFM");


    // Constructors

        //- Construct from a dictionary and a phase pair
        deWildeIPhSATFM
        (
            const dictionary& dict,
            const phasePair& pair,
            const bool registerObject
        );


    //- Destructor
    virtual ~deWildeIPhSATFM();


    // Member Functions

        //- Virtual mass coefficient
        virtual tmp<volScalarField> Cp() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace interPhaseForceModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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
    Foam::dragModels::Tang

(c) Behrad Esgandari 2022
    Department of Particulate Flow Modeling
    Johannes Kepler University, Linz, Austria

Description
    Drag model of Tang et al. for flow past monodisperse static arrays of spheres
    with direct numerical simulations with iterative immersed boundary method.
    range of applicability: 0 < alphaParticles < 0.6, Re <= 1000

    Reference:
    \verbatim
        Tang, Y., E. A. J. F. Peters, J. A. M. Kuipers, S. H. L. Kriebitzsch, 
	and M. A.  van der Hoef. (2015).
        A new drag correlation from fully resolved simulations of 
	flow past monodisperse static arrays of spheres.
        AIChE Journal, 61(2), 688-698.
    \endverbatim

SourceFiles
    Tang.C

\*---------------------------------------------------------------------------*/

#ifndef Tang_H
#define Tang_H

#include "dragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace dragModels
{

/*---------------------------------------------------------------------------*\
                           Class Tang Declaration
\*---------------------------------------------------------------------------*/

class Tang
:
    public dragModel
{
    // Private data

        //- Residual Reynolds Number
        const dimensionedScalar residualRe_;


public:

    //- Runtime type information
    TypeName("Tang");


    // Constructors

        //- Construct from a dictionary and a phase pair
        Tang
        (
            const dictionary& dict,
            const phasePair& pair,
            const bool registerObject
        );


    //- Destructor
    virtual ~Tang();


    // Member Functions

        //- Drag coefficient
        virtual tmp<volScalarField> CdRe() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dragModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

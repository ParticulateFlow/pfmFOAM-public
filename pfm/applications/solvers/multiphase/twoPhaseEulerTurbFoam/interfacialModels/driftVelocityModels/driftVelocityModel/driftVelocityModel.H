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
    Foam::driftVelocityModel

Description

SourceFiles
    driftVelocityModel.C
    newDriftVelocityModel.C

\*---------------------------------------------------------------------------*/

#ifndef driftVelocityModel_H
#define driftVelocityModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class phasePair;

/*---------------------------------------------------------------------------*\
                    Class driftVelocityModel Declaration
\*---------------------------------------------------------------------------*/

class driftVelocityModel
{
protected:

    // Protected data

        //- Phase pair
        const phasePair& pair_;

public:

    //- Runtime type information
    TypeName("driftVelocityModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            driftVelocityModel,
            dictionary,
            (
                const dictionary& dict,
                const phasePair& pair
            ),
            (dict, pair)
        );

    // Static data members

        //- velocitay dimensions
        static const dimensionSet dimU;
        //- Force dimensions
        static const dimensionSet dimF;


    // Constructors

        //- Construct from a dictionary and a phase pair
        driftVelocityModel
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~driftVelocityModel();


    // Selectors

        static autoPtr<driftVelocityModel> New
        (
            const dictionary& dict,
            const phasePair& pair
        );


    // Member Functions

        //- Turbulent drift velocity
        virtual tmp<volVectorField> udrift() const = 0;
        //- Turbulent dispersion force
        virtual tmp<volVectorField> KdUdrift() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

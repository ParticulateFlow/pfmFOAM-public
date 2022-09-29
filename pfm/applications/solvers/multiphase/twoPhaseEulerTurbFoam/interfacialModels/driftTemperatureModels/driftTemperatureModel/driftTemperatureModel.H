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
    Foam::driftTemperatureModel

Description

SourceFiles
    driftTemperatureModel.C
    newDriftTemperatureModel.C

\*---------------------------------------------------------------------------*/

#ifndef driftTemperatureModel_H
#define driftTemperatureModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class phasePair;

/*---------------------------------------------------------------------------*\
                    Class driftTemperatureModel Declaration
\*---------------------------------------------------------------------------*/

class driftTemperatureModel
:
    public regIOobject
{
protected:

    // Protected data

        //- Phase pair
        const phasePair& pair_;
    
        mutable volScalarField heatTransferCorr_;

public:

    //- Runtime type information
    TypeName("driftTemperatureModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            driftTemperatureModel,
            dictionary,
            (
                const dictionary& dict,
                const phasePair& pair,
                const bool registerObject
            ),
            (dict, pair, registerObject)
        );

    // Static data members

        //- velocitay dimensions
        static const dimensionSet dimT;
        //- Force dimensions
        static const dimensionSet dimG;


    // Constructors

        //- Construct from a dictionary and a phase pair
        driftTemperatureModel
        (
            const dictionary& dict,
            const phasePair& pair,
            const bool registerObject
        );


    //- Destructor
    virtual ~driftTemperatureModel();


    // Selectors

        static autoPtr<driftTemperatureModel> New
        (
            const dictionary& dict,
            const phasePair& pair
        );


    // Member Functions

        //- Turbulent drift Temperature
        virtual tmp<volScalarField> Tdrift() const = 0;
        //- Turbulent dispersion force
        virtual tmp<volScalarField> KhTdrift() const;
    
        //- Dummy write for regIOobject
        bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

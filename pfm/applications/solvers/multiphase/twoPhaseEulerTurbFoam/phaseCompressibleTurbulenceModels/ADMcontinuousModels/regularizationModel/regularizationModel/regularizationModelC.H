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
    Foam::ADMcontinuousModels::regularizationModel

Description

SourceFiles
    regularizationModel.C

\*---------------------------------------------------------------------------*/

#ifndef regularizationModelC_H
#define regularizationModelC_H

#include "dictionary.H"
#include "volFields.H"
#include "dimensionedTypes.H"
#include "runTimeSelectionTables.H"
#include "phaseModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace ADMcontinuousModels
{

/*---------------------------------------------------------------------------*\
                           Class regularizationModel Declaration
\*---------------------------------------------------------------------------*/

class regularizationModel
{
    // Private member functions

        //- Disallow default bitwise copy construct
        regularizationModel(const regularizationModel&);

        //- Disallow default bitwise assignment
        void operator=(const regularizationModel&);


protected:

    // Protected data

        const dictionary& dict_;


public:

    //- Runtime type information
    TypeName("regularizationModel");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        regularizationModel,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );


    // Constructors

        //- Construct from components
        regularizationModel(const dictionary& dict);


    // Selectors

        static autoPtr<regularizationModel> New
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~regularizationModel();


    // Member Functions

        virtual tmp<volVectorField> regTerm
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            const volScalarField& k,
            const volVectorField& U,
            const volVectorField& Ustar
        ) const = 0;

        virtual bool read()
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace ADMcontinuousModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

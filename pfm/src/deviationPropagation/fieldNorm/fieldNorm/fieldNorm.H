/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Contributing authors:
    Thomas Lichtenegger
    Copyright (C) 2015- Johannes Kepler University, Linz
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
    Foam::fieldNorm::fieldNorm

Description
    Virtual parent class to compute norms of fields.

\*---------------------------------------------------------------------------*/

#ifndef fieldNorm_H
#define fieldNorm_H

#include "dataBase.H"
#include "HashTable.H"
#include "OFstream.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class fieldNorm Declaration
\*---------------------------------------------------------------------------*/

class fieldNorm
{
protected:

    dataBase& dataBase_;

public:

    friend class dataBase;

    //- Runtime type information
    TypeName("fieldNorm");

    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            fieldNorm,
            dictionary,
            (
                const dictionary& dict,
                dataBase& base
            ),
            (dict,base)
        );


    // Constructors

        //- Construct from components
        fieldNorm
        (
            const dictionary& dict,
            dataBase& base
        );


    // Destructor

        virtual ~fieldNorm();


    // Selector

        static autoPtr<fieldNorm> New
        (
            const dictionary& dict,
            dataBase& base
        );


    // Member functions

    virtual scalar fieldsDistance(const volScalarField &, const volScalarField &, double normalization = -1.0) = 0;

    virtual scalar fieldsDistance(const volVectorField &, const volVectorField &, double normalization = -1.0) = 0;

    virtual scalar fieldsDistance(const volTensorField &, const volTensorField &, double normalization = -1.0) = 0;

    virtual scalar fieldsDistanceConvectiveTerm(const volVectorField &, const volVectorField &, double normalization = -1.0) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

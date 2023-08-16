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
    Foam::fieldNorm::diffNorm

Description
    Class to compute the distance of two fields in terms of their difference.

\*---------------------------------------------------------------------------*/

#ifndef diffNorm_H
#define diffNorm_H

#include "fieldNorm.H"
#include "cellSet.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class diffNorm Declaration
\*---------------------------------------------------------------------------*/

class diffNorm
:
    public fieldNorm
{
protected:

    // Protected data

    dictionary propsDict_;

    bool verbose_;

    word fieldName_;

    scalar domainVolume_;


public:

    //- Runtime type information
    TypeName("diffNorm");

    // Constructors

        //- Construct from components
        diffNorm
        (
            const dictionary& dict,
            dataBase& base,
            word type = "diffNorm"
        );


    // Destructor

        virtual ~diffNorm();


    // Member Functions

    scalar fieldsDistance(const volScalarField &, const volScalarField &, double normalization);

    scalar fieldsDistance(const volVectorField &, const volVectorField &, double normalization);

    scalar fieldsDistance(const volTensorField &, const volTensorField &, double normalization);

    scalar fieldsDistanceConvectiveTerm(const volVectorField &, const volVectorField &, double normalization);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

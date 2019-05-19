/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::ParmentierTestFilter

Description
    Simple top-hat filter used in dynamic LES models with correct, which is
    equivalent to Laplace filter with C = 1/7 in the case of an isotropic
    grid.

    Implemented as a surface integral of the face interpolate of the field.

SourceFiles
    ParmentierTestFilter.C

\*---------------------------------------------------------------------------*/

#ifndef ParmentierTestFilter_H
#define ParmentierTestFilter_H

#include "LESfilter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class ParmentierTestFilter Declaration
\*---------------------------------------------------------------------------*/

class ParmentierTestFilter
:
    public LESfilter

{
    // Private Member Functions

        //- Disallow default bitwise copy construct and assignment
        ParmentierTestFilter(const ParmentierTestFilter&);
        void operator=(const ParmentierTestFilter&);


public:

    //- Runtime type information
    TypeName("ParmentierTest");


    // Constructors

        //- Construct from components
        ParmentierTestFilter(const fvMesh& mesh);

        //- Construct from IOdictionary
        ParmentierTestFilter(const fvMesh& mesh, const dictionary&);


    //- Destructor
    virtual ~ParmentierTestFilter()
    {}


    // Member Functions

        //- Read the LESfilter dictionary
        virtual void read(const dictionary&);


    // Member Operators

        virtual tmp<volScalarField> operator()
        (
            const tmp<volScalarField>&
        ) const;

        virtual tmp<volVectorField> operator()
        (
            const tmp<volVectorField>&
        ) const;

        virtual tmp<volSymmTensorField> operator()
        (
            const tmp<volSymmTensorField>&
        ) const;

        virtual tmp<volTensorField> operator()
        (
            const tmp<volTensorField>&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

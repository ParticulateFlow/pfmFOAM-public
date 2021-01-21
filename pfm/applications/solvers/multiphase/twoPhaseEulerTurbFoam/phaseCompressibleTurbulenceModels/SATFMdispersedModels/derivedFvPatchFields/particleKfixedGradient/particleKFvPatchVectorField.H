/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 
     \\/     M anipulation  | VTT Technical Research Centre of Finland Ltd
-------------------------------------------------------------------------------
License
    This file is a derived work of OpenFOAM.

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
    Foam::particleKFvPatchVectorField
   
Description
    Boundary condition for the TKE of particle phase; use together with
    partial slip boundary condition
 
    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        sigma        | turbulent Schmidt number| no          | 1
        muW          | friction coefficient    | yes         |
    \endtable


    Example usage:
    \verbatim
    roughWall
    {
        type            particleK;
        muW             0.2;
        sigma           1;
        value           $internalField;
    }
    
SourceFiles
    particleKFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef particleKFvPatchVectorField_H
#define particleKFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class particleThetaSoleimaniSchneiderbauerFvPatch Declaration
\*---------------------------------------------------------------------------*/

class particleKFvPatchVectorField
:
    public fixedGradientFvPatchVectorField
{
    // Private data

        //- Wall friction coefficient for non-sliding particles
	    scalar muW_;

        //- Standard deviation (roughness)
        scalar sigma_;

        //- Residual volume fraction
        scalar residualAlpha_;


public:


    //- Runtime type information
    TypeName("particleK");


    // Constructors

        //- Construct from patch and internal field
        particleKFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<Vector<double>, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        particleKFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<Vector<double>, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  particleKFvPatchVectorField
        //  onto a new patch
        particleKFvPatchVectorField
        (
            const particleKFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<Vector<double>, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new particleKFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        particleKFvPatchVectorField
        (
            const particleKFvPatchVectorField&,
            const DimensionedField<Vector<double>, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone(const DimensionedField<Vector<double>, volMesh>& iF) const
        {
            return tmp<fvPatchVectorField>
            (
                new particleKFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchVectorField&,
                const labelList&
            );


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

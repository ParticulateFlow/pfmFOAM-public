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
    Foam::partialSlipFrictionFvPatchVectorField

SourceFiles
    partialSlipFrictionFvPatchVectorField.C
    
Description
    Frictional boundary condition for flat frictional walls.

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        muw          | particle wall friction  | yes         |
        residualAlpha| residual volume fraction | no         | 1e-8
        residualU    | residual velocity       | no          | 1e-8
    \endtable

    Note: If needed,'value' must be defined to be used.

    Example usage:
    \verbatim
    wall
    {
        type            partialSlipFriction;
        muw             0.2;
    }

SourceFiles
    partialSlipFrictionFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef partialSlipFrictionFvPatchVectorField_H
#define partialSlipFrictionFvPatchVectorField_H

#include "fvPatchFields.H"
#include "partialSlipFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class partialSlipFrictionFvPatch Declaration
\*---------------------------------------------------------------------------*/

class partialSlipFrictionFvPatchVectorField
:
    public partialSlipFvPatchVectorField
{
    // Private data

        //- Wall friction coefficient for non-sliding particles
	    scalar muW_;

        //- Residual volume fraction
        scalar residualAlpha_;

        //- Residual velocity
        scalar residualU_;


public:


    //- Runtime type information
    TypeName("partialSlipFriction");


    // Constructors

        //- Construct from patch and internal field
        partialSlipFrictionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<Vector<double>, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        partialSlipFrictionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<Vector<double>, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  partialSlipFrictionFvPatchVectorField onto a new patch
        partialSlipFrictionFvPatchVectorField
        (
            const partialSlipFrictionFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<Vector<double>, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new partialSlipFrictionFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        partialSlipFrictionFvPatchVectorField
        (
            const partialSlipFrictionFvPatchVectorField&,
            const DimensionedField<Vector<double>, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone(const DimensionedField<Vector<double>, volMesh>& iF) const
        {
            return tmp<fvPatchVectorField>
            (
                new partialSlipFrictionFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

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

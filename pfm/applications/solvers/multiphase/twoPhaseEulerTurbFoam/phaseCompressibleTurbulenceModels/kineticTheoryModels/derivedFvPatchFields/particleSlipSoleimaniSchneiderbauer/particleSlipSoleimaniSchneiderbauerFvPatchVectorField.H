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
    Foam::particleSlipSoleimaniSchneiderbauerFvPatchVectorField

SourceFiles
    particleSlipSoleimaniSchneiderbauerFvPatchVectorField.C
    
Description
    Soleimani-Schneiderbauer boundary condition for the velocity of the 
    particles at rough walls.

    References:
    \verbatim
        "Solid boundary condition for collisional gas-solid flows at rough 
        walls"
        Soleimani, A., Pirker, S., Schneiderbauer, S.,
        Powder Technology 281, 2008, pp. 28-22
    \endverbatim

    \verbatim
        "Non-steady state boundary conditions for collisional granular flows
        at flat frictional moving walls"
        Schneiderbauer, S., Schellander, D., Löderer, A., Pirker, S.,
        International Journal of Multiphase Flow 43, 2012, pp. 149-156
    \endverbatim

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        restitutionCoefficient                 | yes         |
        tangentialRestitutionCoeff             | yes         |
        sigma        | Standard deviation (roughness) | yes  |
        residualAlpha| residual volume fraction | no         | 1e-8
        residualTheta| residual granular temperature | no    | 1e-8
        residualU    | residual velocity       | no          | 1e-8
    \endtable

    Note: If needed,'value' must be defined to be used.

    Example usage:
    \verbatim
    roughWall
    {
        type            particleSlipSoleimaniSchneiderbauer;
        restitutionCoefficient 0.9;
        tangentialRestitutionCoeff 0.9;
        sigma           3;
    }

SourceFiles
    particleSlipThetaSoleimaniSchneiderbauerFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef particleSlipSoleimaniSchneiderbauerFvPatchVectorField_H
#define particleSlipSoleimaniSchneiderbauerFvPatchVectorField_H

#include "fvPatchFields.H"
#include "partialSlipFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class particleSlipSoleimaniSchneiderbauerFvPatch Declaration
\*---------------------------------------------------------------------------*/

class particleSlipSoleimaniSchneiderbauerFvPatchVectorField
:
    public partialSlipFvPatchVectorField
{
    // Private data

	    //- Normal restitution coefficent
   	    scalar restitutionCoefficient_;

        //- Tangential restitution coefficient
	    scalar tangentialRestitutionCoeff_;	

        //- Wall friction coefficient for non-sliding particles
	    scalar muF_;

        //- Standard deviation (roughness)
        scalar sigma_;

        //- Residual volume fraction
        scalar residualAlpha_;

        //- Residual granular temperature
        scalar residualTheta_;

        //- Residual velocity
        scalar residualU_;

        // Model parameter calculated from user input
        scalar mu0_;


public:


    //- Runtime type information
    TypeName("particleSlipSoleimaniSchneiderbauer");


    // Constructors

        //- Construct from patch and internal field
        particleSlipSoleimaniSchneiderbauerFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<Vector<double>, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        particleSlipSoleimaniSchneiderbauerFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<Vector<double>, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  particleSlipSoleimaniSchneiderbauerFvPatchVectorField onto a new patch
        particleSlipSoleimaniSchneiderbauerFvPatchVectorField
        (
            const particleSlipSoleimaniSchneiderbauerFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<Vector<double>, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new particleSlipSoleimaniSchneiderbauerFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        particleSlipSoleimaniSchneiderbauerFvPatchVectorField
        (
            const particleSlipSoleimaniSchneiderbauerFvPatchVectorField&,
            const DimensionedField<Vector<double>, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone(const DimensionedField<Vector<double>, volMesh>& iF) const
        {
            return tmp<fvPatchVectorField>
            (
                new particleSlipSoleimaniSchneiderbauerFvPatchVectorField(*this, iF)
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

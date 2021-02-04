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
    Foam::kineticTheoryModels::frictionalStressModels::SchneiderbauerEtAl

Description
    Implementaton of mu(I)-rheology (Schneiderbauer et al., Chem. Eng.
    Science, 80, 2012)
SourceFiles
    SchneiderbauerEtAlFrictionalStress.C

\*---------------------------------------------------------------------------*/

#ifndef SchneiderbauerEtAl_H
#define SchneiderbauerEtAl_H

#include "frictionalStressModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace frictionalStressModels
{

/*---------------------------------------------------------------------------*\
               Class SchneiderbauerEtAl Declaration
\*---------------------------------------------------------------------------*/

class SchneiderbauerEtAl
:
    public frictionalStressModel
{
    // Private data

        dictionary coeffDict_;

        //- Material constant for frictional normal stress
        dimensionedScalar b_;

        //- Material constant for frictional shear stress
        dimensionedScalar muSt_;

        //- Material constant for frictional shear stress
        dimensionedScalar muC_;

        //- inertial number
        dimensionedScalar I0_;
    
        //- alpha_QS
        dimensionedScalar aQSk_;
    
        //- alpha_int
        dimensionedScalar aInt_;
    
        //-spring stiffness
        dimensionedScalar k_;
    
        //- Lower limit for (alphaMax - alpha1)
        dimensionedScalar alphaDeltaMin_;
    
        //- Rolling friction/Cohesion
        dimensionedScalar Rc_;


public:

    //- Runtime type information
    TypeName("SchneiderbauerEtAl");


    // Constructors

        //- Construct from components
        SchneiderbauerEtAl(const dictionary& dict);


    //- Destructor
    virtual ~SchneiderbauerEtAl();


    // Member functions

        virtual tmp<volScalarField> frictionalPressure
        (
            const phaseModel& phase,
            const dimensionedScalar& alphaMinFriction,
            const dimensionedScalar& alphaMax,
            const volScalarField& dp,
            const volScalarField& rho,
            const volSymmTensorField& D
        ) const;

        virtual tmp<volScalarField> frictionalPressurePrime
        (
            const phaseModel& phase,
            const dimensionedScalar& alphaMinFriction,
            const dimensionedScalar& alphaMax,
            const volScalarField& dp,
            const volScalarField& rho,
            const volSymmTensorField& D
        ) const;

        virtual tmp<volScalarField> nu
        (
            const phaseModel& phase,
            const dimensionedScalar& alphaMinFriction,
            const dimensionedScalar& alphaMax,
            const volScalarField& pf,
            const volScalarField& dp,
            const volSymmTensorField& D
        ) const;

        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace frictionalStressModels
} // End namespace kineticTheoryModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

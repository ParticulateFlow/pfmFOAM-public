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
    Foam::kineticTheoryModels::frictionalStressModels::ChialvoEtAl

 Description
     Combination of Schneiderbauer et al. and Chialvo et al. frictional model

     Reference:
     \verbatim
        Chialvo, S., J. Sun, S. Sundaresan. Phys. Rev. E, 2012, 85, 021305 (2012).
        Schneiderbauer, S., A. Aigner, S. Pirker. Chem. Eng. Sci., 2012, 80, 279–292.
     \endverbatim
SourceFiles
    ChialvoEtAlFrictionalStress.C

\*---------------------------------------------------------------------------*/

#ifndef ChialvoEtAl_H
#define ChialvoEtAl_H

#include "frictionalStressModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace frictionalStressModels
{

/*---------------------------------------------------------------------------*\
               Class ChialvoEtAl Declaration
\*---------------------------------------------------------------------------*/

class ChialvoEtAl
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
    TypeName("ChialvoEtAl");


    // Constructors

        //- Construct from components
        ChialvoEtAl(const dictionary& dict);


    //- Destructor
    virtual ~ChialvoEtAl();


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

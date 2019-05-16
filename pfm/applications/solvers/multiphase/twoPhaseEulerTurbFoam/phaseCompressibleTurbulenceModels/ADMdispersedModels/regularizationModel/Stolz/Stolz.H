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
    Foam::ADMdispersedModels::regularizationModels::Stolz
 
(c) Simon Schneiderbauer 2019
    Christian Doppler Laboratory for Multiscale Modeling of Multiphase Processes
    Johannes Kepler University, Linz, Austria

Description
    Regularization model for Approximate deconvolution model

    Reference:
    \verbatim
        Schneiderbauer, S., & Saeedipour, M. (2019). Numerical simulation of turbulent gas-solid flow using an approximate deconvolution model. International Journal of Multiphase Flow, 114, 287–302.
    \endverbatim
 
SourceFiles
    Stolz.C

\*---------------------------------------------------------------------------*/

#ifndef Stolz_H
#define Stolz_H

#include "regularizationModel.H"
// h-files for LES filters
#include "LESfilter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace ADMdispersedModels
{
namespace regularizationModels
{

/*---------------------------------------------------------------------------*\
                           Class Stolz Declaration
\*---------------------------------------------------------------------------*/

class Stolz
:
    public regularizationModel
{
    dictionary coeffDict_;
    //- smagorinsky constant according to Schneiderbauer (2017)
    dimensionedScalar smagConst_;
    
    //- length scale constant for turbulent time scale
    dimensionedScalar lengthConst_;
    
protected:
    autoPtr<LESfilter> filterPtr_;
    LESfilter& filter_;
    
public:

    //- Runtime type information
    TypeName("Stolz");

    // Constructors

        Stolz
        (
            const dictionary& dict,
            const volScalarField& alpha
         );


    //- Destructor
    virtual ~Stolz();


    // Member functions

        virtual tmp<volVectorField> regTerm
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            const volScalarField& k,
            const volVectorField& U,
            const volVectorField& Ustar
        ) const;
    
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace regularizationModels
} // End namespace ADMdispersedModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

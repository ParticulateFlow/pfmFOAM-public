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
    Foam::RASModels::ADMcontinuousModel

(c) Simon Schneiderbauer 2019
    Christian Doppler Laboratory for Multiscale Modeling of Multiphase Processes
    Johannes Kepler University, Linz, Austria

Description
    Approximate deconvolution model for  continuous phase

    Reference:
    \verbatim
        Schneiderbauer, S., & Saeedipour, M. (2019). Numerical simulation of turbulent gas-solid flow using an approximate deconvolution model. International Journal of Multiphase Flow, 114, 287–302.
        Schneiderbauer, S., & Saeedipour, M. (2018). Approximate deconvolution model for the simulation of turbulent gas-solid flows: An a-priori analysis. Physics of Fluids, 30(2), 023301.
    \endverbatim

    There are no default model coefficients.

SourceFiles
    ADMcontinuousModel.C

\*---------------------------------------------------------------------------*/

#ifndef ADMcontinuousModel_H
#define ADMcontinuousModel_H

#include "RASModel.H"
#include "eddyViscosity.H"
#include "phaseCompressibleTurbulenceModel.H"
#include "EddyDiffusivity.H"
#include "phaseModel.H"
#include "dragModel.H"
#include "regularizationModelC.H"
// h-files for LES filters
#include "LESfilter.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class ADMcontinuousModel Declaration
\*---------------------------------------------------------------------------*/

class ADMcontinuousModel
:
    public eddyViscosity
    <
        RASModel<EddyDiffusivity<phaseCompressibleTurbulenceModel>>
    >
{
    // Private data

        // Input Fields

            const phaseModel& phase_;


        // Sub-models
   
            //- Run-time selected regularization model
            autoPtr<ADMcontinuousModels::regularizationModel> regularizationModel_;

        // ADM dispersed coefficients

            //- Maximum packing phase-fraction
            dimensionedScalar alphaMax_;

            //- Residual phase fraction
            dimensionedScalar residualAlpha_;

            //- Maximum turbulent kinetic energy
            dimensionedScalar maxK_;
    
            //- deconvolution order
            dimensionedScalar deconOrder_;
    
            //- Maximum turbulent viscosity
            dimensionedScalar maxNut_;
    
            //- Smagorinsky constant to compute nut_
            dimensionedScalar Cmu_;


        // ADM Fields Fields

            //- Reynolds stress tensor
            volTensorField R2ADM_;
    
            //- deconvolution of volume fraction
            volScalarField alpha2star_;
    
            //- deconvolution of solids velocity
            volVectorField U2star_;
    
            //- Turbulent kinetic energy
            volScalarField k_;


    // Private Member Functions

        void correctNut()
        {}

        //- Disallow default bitwise copy construct
        ADMcontinuousModel(const ADMcontinuousModel&);

        //- Disallow default bitwise assignment
        void operator=(const ADMcontinuousModel&);

protected:
    autoPtr<LESfilter> filterPtr_;
    LESfilter& filter_;
    
    //- bound k
    void boundNormalStress(volTensorField& R) const;

public:

    //- Runtime type information
    TypeName("ADMcontinuous");


    // Constructors

        //- Construct from components
        ADMcontinuousModel
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const phaseModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~ADMcontinuousModel();


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the effective viscosity
        virtual tmp<volScalarField> nuEff() const
        {
            return this->nut() + this->nu()*alpha_;
        }

        //- Return the effective viscosity on patch
        virtual tmp<scalarField> nuEff(const label patchi) const
        {
            return this->nut(patchi) + this->nu(patchi)*alpha_;
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const;

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const;

        //- Return the Reynolds stress tensor
        virtual tmp<volSymmTensorField> R() const;

        //- Return the effective stress tensor
        virtual tmp<volSymmTensorField> devRhoReff() const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevRhoReff(volVectorField& U) const;

        //- Solve the ADM equations and correct the viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    Foam::RASModels::ADMdispersedModel

(c) Simon Schneiderbauer 2019
    Christian Doppler Laboratory for Multiscale Modeling of Multiphase Processes
    Johannes Kepler University, Linz, Austria

Description
    Approximate deconvolution model for continuous phase

    Reference:
    \verbatim
        Schneiderbauer, S., & Saeedipour, M. (2019). Numerical simulation of turbulent gas-solid flow using an approximate deconvolution model. International Journal of Multiphase Flow, 114, 287–302.
        Schneiderbauer, S., & Saeedipour, M. (2018). Approximate deconvolution model for the simulation of turbulent gas-solid flows: An a-priori analysis. Physics of Fluids, 30(2), 023301.
    \endverbatim

    There are no default model coefficients.

SourceFiles
    ADMdispersedModel.C

\*---------------------------------------------------------------------------*/

#ifndef ADMdispersedModel_H
#define ADMdispersedModel_H

#include "RASModel.H"
#include "eddyViscosity.H"
#include "phaseCompressibleTurbulenceModel.H"
#include "EddyDiffusivity.H"
#include "phaseModel.H"
#include "dragModel.H"
#include "frictionalStressModelADM.H"
#include "regularizationModel.H"
// h-files for LES filters
#include "LESfilter.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class ADMdispersedModel Declaration
\*---------------------------------------------------------------------------*/

class ADMdispersedModel
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
            autoPtr<ADMdispersedModels::regularizationModel> regularizationModel_;

            //- Run-time selected frictional stress model
            autoPtr<ADMdispersedModels::frictionalStressModel>
                frictionalStressModel_;    

        // ADM dispersed coefficients

            //- Maximum packing phase-fraction
            dimensionedScalar alphaMax_;

            //- Min value for which the frictional stresses are zero
            dimensionedScalar alphaMinFriction_;

            //- Residual phase fraction
            dimensionedScalar residualAlpha_;

            //- Maximum turbulent kinetic energy
            dimensionedScalar maxK_;
    
            //- deconvolution order
            dimensionedScalar deconOrder_;
    
            //- Maximum turbulent viscosity
            dimensionedScalar maxNut_;    


        // ADM Fields Fields

            //- Reynolds stress tensor
            volTensorField R1ADM_;
    
            //- velocity gradient
            volSymmTensorField Dstar_;
    
            //- deconvolution of volume fraction
            volScalarField alpha1star_;

            //- deconvolution variance of solid volume fraction
            volScalarField alphaP2Mean_;
    
            //- deconvolution of solids velocity
            volVectorField U1star_;

            //- The frictional viscosity
            volScalarField nuFric_;
    
            //- The frictional pressure
            volScalarField pf_;
    
            //- Turbulent kinetic energy
            volScalarField k_;


    // Private Member Functions

        void correctNut()
        {}

        //- Disallow default bitwise copy construct
        ADMdispersedModel(const ADMdispersedModel&);

        //- Disallow default bitwise assignment
        void operator=(const ADMdispersedModel&);

protected:
    autoPtr<LESfilter> filterPtr_;
    LESfilter& filter_;
    
    //- bound k
    void boundNormalStress(volTensorField& R) const;

public:

    //- Runtime type information
    TypeName("ADMdispersed");


    // Constructors

        //- Construct from components
        ADMdispersedModel
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
    virtual ~ADMdispersedModel();


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the effective viscosity
        virtual tmp<volScalarField> nuEff() const
        {
            return this->nut();
        }

        //- Return the effective viscosity on patch
        virtual tmp<scalarField> nuEff(const label patchi) const
        {
            return this->nut(patchi);
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const;

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const;

        //- Return the Reynolds stress tensor
        virtual tmp<volSymmTensorField> R() const;

        //- Return the phase-pressure'
        // (derivative of phase-pressure w.r.t. phase-fraction)
        virtual tmp<volScalarField> pPrime() const;

        //- Return the face-phase-pressure'
        // (derivative of phase-pressure w.r.t. phase-fraction)
        virtual tmp<surfaceScalarField> pPrimef() const;
    
        //- Return the phase-pressure'
        virtual tmp<volVectorField> divStress() const;

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

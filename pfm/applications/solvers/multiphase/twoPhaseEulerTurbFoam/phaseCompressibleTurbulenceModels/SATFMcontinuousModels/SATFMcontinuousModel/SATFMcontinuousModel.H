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
    Foam::RASModels::SATFMcontinuousModel
 
(c) Simon Schneiderbauer 2019
    Christian Doppler Laboratory for Multiscale Modeling of Multiphase Processes
    Johannes Kepler University, Linz, Austria 

Description
    SA-TFM Model

    Reference:
    \verbatim
        Schneiderbauer, S. (2017). A spatially-averaged two-fluid model for dense large-scale gas-solid flows. AIChE Journal, 63(8), 3544–3562.
        Schneiderbauer, S. (2018). Validation study on spatially averaged two-fluid model for gas-solid flows: I. A-priori analysis of wall bounded flows. AIChE Journal, 64(5), 1591–1605.
        Rauchenzauner, S. and Schneiderbauer, S. (2020). A Dynamic Anisotropic Spatially-Averaged Two-Fluid Model for Moderately Dense Gas-Particle Flows. IJMF, 126, 103237.
        Rauchenzauner, S. and Schneiderbauer, S. (2022). A dynamic multiphase turbulence model for coarse-grid simulations of fluidized gas-particle suspensions. Chem. Eng. Sci., 1247, 117104.
    \endverbatim

SourceFiles
    SATFMcontinuousModel.C

\*---------------------------------------------------------------------------*/

#ifndef SATFMcontinuousModel_H
#define SATFMcontinuousModel_H

#include "RASModel.H"
#include "eddyViscosity.H"
#include "phaseCompressibleTurbulenceModel.H"
#include "EddyDiffusivity.H"
#include "phaseModel.H"
#include "dragModel.H"
#include "driftVelocityModel.H"
#include "frictionalStressModel.H"
// h-files for LES filters
#include "LESfilter.H"
#include "simpleFilter.H"
#include "simpleFilterADM.H"
#include "laplaceFilter.H"
#include "wallFvPatch.H"
#include "cyclicAMIFvsPatchField.H"
#include "cyclicACMIFvsPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class SATFMcontinuousModel Declaration
\*---------------------------------------------------------------------------*/

class SATFMcontinuousModel
:
    public eddyViscosity
    <
        RASModel<EddyDiffusivity<phaseCompressibleTurbulenceModel>>
    >
{
    // Private data

        // Input Fields

            const phaseModel& phase_;


        // SA-TFM Model coefficients

            //- Use equilibrium approximation: generation == dissipation
            Switch equilibriumK_;
            Switch equilibriumPhiP2_;
            Switch equilibriumH_;

            //- Use dynamic adjustment of several model parameters and
            //  correlation coefficients
            Switch dynamicAdjustment_;
    
            //- switch between isotropic and anisotropic treatment of
            //  shear component of mesoscale stress tensor
            Switch anIsoTropicNut_;
    
            //- Maximum packing phase-fraction
            dimensionedScalar alphaMax_;    

            //- Residual phase fraction
            dimensionedScalar residualAlpha_;

            //- Maximum turbulent viscosity
            dimensionedScalar maxNut_;
    
            //- correlation coefficient between solid volume fraction and cont. velocity
            dimensionedScalar xiPhiContScalar_;

            //- correlation coefficient between solid volume fraction and cont. temperature
            dimensionedScalar xiTPhiGScalar_;

            //- "Smagorinsky constant"
            dimensionedScalar CmuScalar_;
    
            //- "Smagorinsky constant" near walls
            dimensionedScalar CmuWScalar_;
    
            //- "dissipation constant"
            dimensionedScalar CepsScalar_;
    
            //- "pressure dilation constant"
            dimensionedScalar CpScalar_;

            //- "pressure transport constant"
            dimensionedScalar Cdp_;
    
            //- turbulent schmidt number
            dimensionedScalar sigma_;
    
            //- gravity vector
            dimensionedVector gN_;
    
            //- maximum turbulent kinetic energy
            dimensionedScalar maxK_;
    
        // SA-TFM fields

            //- normal components of Reynoldsstress tensor
            //  (turbulent kinetic energies)
            volVectorField k_;

            //- turbulent internal energy
            volScalarField H_;
    
            //- correlation coefficient between solid volume fraction and cont. velocity
            volVectorField xiPhiG_;

            //- pressure dilatation
            //volScalarField pdil_;

            //- correlation coefficient between solid volume fraction and cont. temperature
            volScalarField xiTPhiG_;
    
            //- auto correlation coefficient of solids velocity
            volTensorField xiUU_;
    
            //- smagorinsky constant
            volScalarField Cmu_;
    
            //- turbulent dissipation constant
            volScalarField Ceps_;
    
            //- pressure dilation constant
            volScalarField Cp_;
    
            //- grid size
            volScalarField deltaF_;
    
            //- mixing length
            volScalarField lm_;
    
            //- Reynolds stress tensor
            volTensorField R2_;
    
            //- anisotropic viscosity
            volTensorField nutA_;
    

    // Private Member Functions

        void correctNut()
        {}

        //- Disallow default bitwise copy construct
        SATFMcontinuousModel(const SATFMcontinuousModel&);

        //- Disallow default bitwise assignment
        void operator=(const SATFMcontinuousModel&);

protected:
    autoPtr<LESfilter> filterPtr_;
    LESfilter& filter_;
    
    //- bound k
    void boundNormalStress(volVectorField& R) const;
    void boundxiPhiG(volVectorField& xi) const;
    void boundCorrTensor(volTensorField& R) const;
    void boundGradU(volTensorField& S) const;
    void boundStress(volTensorField& R) const;
    void boundgradh(volVectorField& gradh) const;
public:

    //- Runtime type information
    TypeName("SATFMcontinuous");


    // Constructors

        //- Construct from components
        SATFMcontinuousModel
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
    virtual ~SATFMcontinuousModel();


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
    
        //- Return the phase-pressure'
        virtual tmp<volScalarField> normalStress() const;

        //- Solve the SA-TFM equations and correct the viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2018 OpenFOAM Foundation
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
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/
 
  (c) Simon Schneiderbauer 2019
    Christian Doppler Laboratory for Multiscale Modeling of Multiphase Processes
    Johannes Kepler University, Linz, Austria

\*---------------------------------------------------------------------------*/

#include "driftVelocityADMreg.H"
#include "phasePair.H"
#include "PhaseCompressibleTurbulenceModel.H"
#include "addToRunTimeSelectionTable.H"

#include "dragModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace driftVelocityModels
{
    defineTypeNameAndDebug(driftVelocityADMreg, 0);
    addToRunTimeSelectionTable
    (
        driftVelocityModel,
        driftVelocityADMreg,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::driftVelocityModels::driftVelocityADMreg::driftVelocityADMreg
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject
)
:
    driftVelocityModel(dict, pair, registerObject),
    residualAlpha_
    (
        "residualAlpha",
        dimless,
        dict.lookupOrDefault<scalar>
        (
            "residualAlpha",
            pair_.dispersed().residualAlpha().value()
        )
    ),

    alphaMax_("alphaMax", dimless, dict),

    filterPtr_(LESfilter::New(pair.dispersed().mesh(), dict)),
    filter_(filterPtr_())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::driftVelocityModels::driftVelocityADMreg::~driftVelocityADMreg()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField>
Foam::driftVelocityModels::driftVelocityADMreg::udrift() const
{
    const fvMesh& mesh_(pair_.phase1().mesh());
    const volScalarField& alpha1star(mesh_.lookupObject<volScalarField>
                               ("alpha1star"));
    
    const volVectorField& U2star(mesh_.lookupObject<volVectorField>
                           ("U2star"));
    const volScalarField& alphaP2Mean_(mesh_.lookupObject<volScalarField>
                                       ("alphaP2Mean"));

    volScalarField alpha1f = filter_(alpha1star);
    alpha1f.max(residualAlpha_.value());
    alpha1f.min(alphaMax_.value());
    volScalarField alpha2f = scalar(1.0) - alpha1f;
    
    volVectorField uSlipV(pair_.continuous().U() - pair_.dispersed().U());
    volScalarField uSlip(mag(uSlipV));
    uSlip.max(SMALL);
    
    return pos(pair_.dispersed() - residualAlpha_)*
           (
                filter_(alpha1star*U2star)/alpha1f
              - filter_((scalar(1.0) - alpha1star)*U2star)/alpha2f
              - sqrt(2.0*alphaP2Mean_*pair_.continuous().turbulence().k())*uSlipV
                    /(uSlip
                        * (
                             max(pair_.dispersed(),residualAlpha_)
                           * (scalar(1.0) - pair_.dispersed())
                          )
                      )
           );
}


// ************************************************************************* //

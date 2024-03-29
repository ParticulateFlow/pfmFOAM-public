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
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
 
 (c) Simon Schneiderbauer 2019
    Christian Doppler Laboratory for Multiscale Modeling of Multiphase Processes
    Johannes Kepler University, Linz, Austria

\*---------------------------------------------------------------------------*/

#include "deWildeIPhADM.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace interPhaseForceModels
{
    defineTypeNameAndDebug(deWildeIPhADM, 0);
    addToRunTimeSelectionTable
    (
        interPhaseForceModel,
        deWildeIPhADM,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interPhaseForceModels::deWildeIPhADM::deWildeIPhADM
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject
)
:
    interPhaseForceModel(dict, pair, registerObject),
    residualAlpha_
    (
        "residualAlpha",
        dimless,
        dict.lookupOrDefault<scalar>
        (
            "residualAlpha",
            pair_.dispersed().residualAlpha().value()
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interPhaseForceModels::deWildeIPhADM::~deWildeIPhADM()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::interPhaseForceModels::deWildeIPhADM::Cp() const
{
    // get alphaP2Mean from Turbulence Model
    const fvMesh& mesh_(pair_.phase1().mesh());
    const volScalarField& alphaP2Mean_(mesh_.lookupObject<volScalarField>
                               ("alphaP2Mean"));
    volScalarField alpha1
    (
        max(pair_.dispersed(), residualAlpha_)
    );
    volScalarField alpha2
    (
        scalar(1.0) - alpha1
    );
    
    // limit alphaP2Mean to prevent unphysical values of Vm0
    volScalarField alphaP2Mean = min(0.9*sqr(alpha1),alphaP2Mean_);
    
    volScalarField rho1
    (
        pair_.dispersed().rho()
    );
    volScalarField rho2
    (
        pair_.continuous().rho()
    );

    volScalarField rho
    (
        alpha1*rho1 + alpha2*rho2
    );
    volScalarField oneMinMg
    (
        alpha2*(scalar(1.0) - rho2/rho)
    );
    
    oneMinMg.max(0);

    volScalarField Cp
    (
        min(alphaP2Mean*rho1/rho,oneMinMg)
    );

    return
        pos(pair_.dispersed() - residualAlpha_)*Cp;
}


// ************************************************************************* //

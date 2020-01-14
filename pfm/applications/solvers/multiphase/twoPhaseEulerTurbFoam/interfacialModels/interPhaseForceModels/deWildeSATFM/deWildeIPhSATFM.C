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

#include "deWildeIPhSATFM.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace interPhaseForceModels
{
    defineTypeNameAndDebug(deWildeIPhSATFM, 0);
    addToRunTimeSelectionTable
    (
        interPhaseForceModel,
        deWildeIPhSATFM,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interPhaseForceModels::deWildeIPhSATFM::deWildeIPhSATFM
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

Foam::interPhaseForceModels::deWildeIPhSATFM::~deWildeIPhSATFM()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::interPhaseForceModels::deWildeIPhSATFM::Cp() const
{
    // get alphaP2Mean from Turbulence Model
    const fvMesh& mesh(pair_.phase1().mesh());
    const volScalarField& alphaP2Mean1_(mesh.lookupObject<volScalarField>
                                        ("alphaP2Mean." + pair_.dispersed().name()));

    const volScalarField& alphaP2Mean2_(mesh.lookupObject<volScalarField>
                                        ("alphaP2Mean." + pair_.continuous().name()));
    
    volScalarField alphaP2Mean = max(alphaP2Mean1_,alphaP2Mean2_);

    volScalarField alpha1
    (
        max(pair_.dispersed(), residualAlpha_)
    );
    volScalarField alpha2
    (
        scalar(1.0) - alpha1
    );
    
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
        alpha2*(scalar(1.0) - alpha2*rho2/rho)
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

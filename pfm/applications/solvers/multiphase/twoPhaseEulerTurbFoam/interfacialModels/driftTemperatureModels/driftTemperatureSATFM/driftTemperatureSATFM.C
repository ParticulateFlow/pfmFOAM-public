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
 
  (c) 2022
    Christian Doppler Laboratory for Multiscale Modeling of Multiphase Processes
    Johannes Kepler University, Linz, Austria

\*---------------------------------------------------------------------------*/

#include "driftTemperatureSATFM.H"
#include "phasePair.H"
#include "PhaseCompressibleTurbulenceModel.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace driftTemperatureModels
{
    defineTypeNameAndDebug(driftTemperatureSATFM, 0);
    addToRunTimeSelectionTable
    (
        driftTemperatureModel,
        driftTemperatureSATFM,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::driftTemperatureModels::driftTemperatureSATFM::driftTemperatureSATFM
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject
)
:
    driftTemperatureModel(dict, pair, registerObject),
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

Foam::driftTemperatureModels::driftTemperatureSATFM::~driftTemperatureSATFM()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::driftTemperatureModels::driftTemperatureSATFM::Tdrift() const
{
    const fvMesh& mesh = pair_.phase1().mesh();
    const volScalarField& alphaP2Mean1_ = mesh.lookupObject<volScalarField>("alphaP2Mean." + pair_.dispersed().name());

    const volScalarField& xiTPhiG_= mesh.lookupObject<volScalarField>("xiTPhiG");

    const volScalarField& HC_ = mesh.lookupObject<volScalarField>("H." + pair_.continuous().name());

    volScalarField alphaP2Mean = alphaP2Mean1_;

    volScalarField alpha1 = max(pair_.dispersed(), residualAlpha_);
    
    volScalarField alphaP2MeanN = sqrt(alphaP2Mean)
                                /(alpha1*(scalar(1.0)-alpha1));

    const rhoThermo& thermo2 = pair_.continuous().thermo();
    //cv or cp
    volScalarField Cpv2("Cpv2", thermo2.Cpv());

    return //pos(pair_.dispersed() - residualAlpha_)*
          alphaP2MeanN
         * sqrt((HC_))
         * xiTPhiG_
         /Cpv2;

}

// ************************************************************************* //

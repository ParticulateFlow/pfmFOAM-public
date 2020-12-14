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

#include "driftVelocitySATFMreg.H"
#include "phasePair.H"
#include "PhaseCompressibleTurbulenceModel.H"
#include "addToRunTimeSelectionTable.H"

#include "dragModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace driftVelocityModels
{
    defineTypeNameAndDebug(driftVelocitySATFMreg, 0);
    addToRunTimeSelectionTable
    (
        driftVelocityModel,
        driftVelocitySATFMreg,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::driftVelocityModels::driftVelocitySATFMreg::driftVelocitySATFMreg
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
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::driftVelocityModels::driftVelocitySATFMreg::~driftVelocitySATFMreg()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField>
Foam::driftVelocityModels::driftVelocitySATFMreg::udrift() const
{
    const fvMesh& mesh(pair_.phase1().mesh());
    const volScalarField& alphaP2Mean(mesh.lookupObject<volScalarField>
                               ("alphaP2Mean." + pair_.dispersed().name()));

    const volVectorField& xiPhiG_(mesh.lookupObject<volVectorField>
                               ("xiPhiG"));
    const volVectorField& kC_(mesh.lookupObject<volVectorField>
                                 ("k." + pair_.continuous().name()));
    dimensionedVector eSum
    (
        "eX",
        dimensionSet(0, 0, 0, 0, 0, 0, 0),
        vector(1,1,1)
    );

    volScalarField alpha1 = max(pair_.dispersed(), residualAlpha_);
    
    volVectorField kSqrt =  xiPhiG_*(kC_&eSum);
    
    volScalarField alphaP2MeanN = sqrt(2.0 * alphaP2Mean)
                                /(alpha1*(scalar(1.0)-alpha1));
    
    return pos(pair_.dispersed() - residualAlpha_)
         * alphaP2MeanN
         * kSqrt;

}


// ************************************************************************* //

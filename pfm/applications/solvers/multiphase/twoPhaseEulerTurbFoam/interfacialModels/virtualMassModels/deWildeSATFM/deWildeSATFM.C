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

#include "deWildeSATFM.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace virtualMassModels
{
    defineTypeNameAndDebug(deWildeSATFM, 0);
    addToRunTimeSelectionTable
    (
        virtualMassModel,
        deWildeSATFM,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::virtualMassModels::deWildeSATFM::deWildeSATFM
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject
)
:
    virtualMassModel(dict, pair, registerObject),
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

Foam::virtualMassModels::deWildeSATFM::~deWildeSATFM()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::virtualMassModels::deWildeSATFM::Cvm() const
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
        pair_.dispersed().rho()
    );

    volScalarField rho
    (
        alpha1*rho1 + alpha2*rho2
    );
    volScalarField Vm0
    (
        alphaP2Mean/(max(alpha1*alpha2-alphaP2Mean,sqr(residualAlpha_)))
       *((alpha1*alpha2*rho1*rho2)/(rho*rho))
    );
    // Limit virtual mass coefficient
    Vm0.min(100.0);
    Vm0.max(0.0);
    
    return
        Vm0*rho/(alpha1*rho2);
}


// ************************************************************************* //

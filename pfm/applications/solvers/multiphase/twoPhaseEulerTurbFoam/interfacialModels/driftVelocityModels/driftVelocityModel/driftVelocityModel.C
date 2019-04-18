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

\*---------------------------------------------------------------------------*/

#include "driftVelocityModel.H"
#include "phasePair.H"
#include "dragModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(driftVelocityModel, 0);
    defineRunTimeSelectionTable(driftVelocityModel, dictionary);
}

const Foam::dimensionSet Foam::driftVelocityModel::dimU(0, 1, -1, 0, 0);
const Foam::dimensionSet Foam::driftVelocityModel::dimF(1, -2, -2, 0, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::driftVelocityModel::driftVelocityModel
(
    const dictionary& dict,
    const phasePair& pair
)
:
    pair_(pair)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::driftVelocityModel::~driftVelocityModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField>
Foam::driftVelocityModel::KdUdrift() const
{
    const fvMesh& mesh(pair_.phase1().mesh());
    const dragModel&
    drag
    (
        mesh.lookupObject<dragModel>
        (
            IOobject::groupName(dragModel::typeName, pair_.name())
        )
    );
    volVectorField ud(udrift());

    volScalarField Re(pair_.Re());
    volScalarField uslip = Re*pair_.continuous().nu()/(pair_.dispersed().d());
     
    // limit turbulent dispersion force according to
    // Parmentier et al., AIChE J., 2012
    
    volScalarField magUd = mag(ud);
    magUd.max(SMALL);
    ud *= min(0.9*uslip,mag(ud))/magUd;
    volScalarField dragCorr = mag(ud)/uslip;
    /*
    Info << "Drift Velocity ADM:" << nl
         << "max. drag correction: " << max(dragCorr).value() << nl
         << "max. drift velocity:  " << max(mag(ud)).value()  << endl;
    */
    // multiply drift velocity by drag coefficient
    return
        0.75
        *pair_.dispersed()
        *drag.CdRe()
        *pair_.continuous().nu()
        *pair_.continuous().rho()
        *ud
        /sqr(pair_.dispersed().d());
}

// ************************************************************************* //

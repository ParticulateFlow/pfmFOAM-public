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
#include "zeroGradientFvPatchField.H"

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
    const phasePair& pair,
    const bool registerObject
)
:
    regIOobject
    (
        IOobject
        (
            IOobject::groupName(typeName, pair.name()),
            pair.phase1().mesh().time().timeName(),
            pair.phase1().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            registerObject
        )
    ),

    pair_(pair),

    dragCorr_
    (
        IOobject
        (
            "dragCorr",
            pair.dispersed().time().timeName(),
            pair.dispersed().mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pair.dispersed().mesh(),
        dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.0),
        // Set Boundary condition
        zeroGradientFvPatchField<scalar>::typeName
        
    )
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
    dimensionedScalar uSmall("uSmall", dimensionSet(0, 1, -1, 0, 0, 0, 0), 1.0e-6);
    
    dimensionedVector eX
    (
        "eX",
        dimensionSet(0, 0, 0, 0, 0, 0, 0),
        vector(1,0,0)
    );
    dimensionedVector eY
    (
        "eY",
        dimensionSet(0, 0, 0, 0, 0, 0, 0),
        vector(0,1,0)
    );
    dimensionedVector eZ
    (
        "eZ",
        dimensionSet(0, 0, 0, 0, 0, 0, 0),
        vector(0,0,1)
    );
    
    volVectorField ud(udrift());
    
    volVectorField uSlipV(pair_.continuous().U() - pair_.dispersed().U());
    volScalarField uSlip(mag(uSlipV));
    uSlip.max(SMALL);
    
    dragCorr_ = -((ud&uSlipV)/sqr(uSlip));
    dragCorr_.min(0.999);
    dragCorr_.max(-0.999);
    
    // limit turbulent dispersion force according to
    // Parmentier et al., AIChE J., 2012
    volScalarField magUd = mag(ud);
    magUd.max(SMALL);
    ud *= min(mag(uSlipV&ud)/uSlip,0.999*uSlip)/magUd;
    /*
    ud =  ((ud&eX)*min(0.99*mag(uSlipV&eX)/(mag(ud&eX)+uSmall),1.0))*eX
        + ((ud&eY)*min(0.99*mag(uSlipV&eY)/(mag(ud&eY)+uSmall),1.0))*eY
        + ((ud&eZ)*min(0.99*mag(uSlipV&eZ)/(mag(ud&eZ)+uSmall),1.0))*eZ;
    */
    // multiply drift velocity by drag coefficient
    return
         drag.K()
        *ud; 
}

bool Foam::driftVelocityModel::writeData(Ostream& os) const
{
    return os.good();
}

// ************************************************************************* //

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

#include "driftTemperatureModel.H"
#include "phasePair.H"
#include "heatTransferModel.H"
#include "zeroGradientFvPatchField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(driftTemperatureModel, 0);
    defineRunTimeSelectionTable(driftTemperatureModel, dictionary);
}

const Foam::dimensionSet Foam::driftTemperatureModel::dimT(0, 0, 0, 1, 0);
const Foam::dimensionSet Foam::driftTemperatureModel::dimG(1, -1, -3, 0, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::driftTemperatureModel::driftTemperatureModel
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

    heatTransferCorr_
    (
        IOobject
        (
            "heatTransferCorr",
            pair.dispersed().time().timeName(),
            pair.dispersed().mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pair.dispersed().mesh(),
        dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.0),
        // Set Boundary condition
        zeroGradientFvPatchField<scalar>::typeName
        
    ),

    driftTemp_
    (
        IOobject
        (
            "driftTemp",
            pair.dispersed().time().timeName(),
            pair.dispersed().mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pair.dispersed().mesh(),
        dimensionedScalar("zero", dimensionSet(0, 0, 0, 1, 0), 0.0),
        // Set Boundary condition
        zeroGradientFvPatchField<scalar>::typeName

    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::driftTemperatureModel::~driftTemperatureModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::driftTemperatureModel::KhTdrift() const
{

    const fvMesh& mesh(pair_.phase1().mesh());
    const heatTransferModel&
    heatTransfer
    (
        mesh.lookupObject<heatTransferModel>
        (
            IOobject::groupName(heatTransferModel::typeName, pair_.name())
        )
    );


    dimensionedScalar TSmall("TSmall", dimensionSet(0, 0, 0, 1, 0, 0, 0), 1.0e-6);
    
    volScalarField Td(Tdrift());

    volScalarField alpha1 = max(pair_.dispersed(), SMALL);

    const rhoThermo& thermo2 = pair_.continuous().thermo();
    //cv or cp and enthalpie h or internal energy e
    volScalarField Cpv2("Cpv2", thermo2.Cpv());
    const volScalarField& heC = thermo2.he();

    const rhoThermo& thermo1 = pair_.dispersed().thermo();
    //cv or cp and enthalpie h or internal energy e of other phase
    volScalarField Cpv1("Cpv1", thermo1.Cpv());
    const volScalarField& heD = thermo1.he();

    //continuous phase temperature
    volScalarField TC = (heC/Cpv2);
    volScalarField TD = (heD/Cpv1);

    volScalarField TSlip(TC - TD);
    volScalarField TSlipM(mag(TSlip));
    TSlipM.max(TSmall);


    //Td = (Td*min((0.99999999999*mag(TSlip)/(mag(Td)+TSmall)),1.0));
    
    heatTransferCorr_ = -(TSlip)*Td/(TSlipM*TSlipM);
    //heatTransferCorr_ = - pos(mag(TSlip) - mag(Td))*Td/(TSlip +TSmall) + neg(mag(TSlip) - mag(Td))*0.99;


    Td = (Td*min((0.99999999999*mag(TSlip)/(mag(Td)+TSmall)),1.0));

    heatTransferCorr_.min(0.99);
    heatTransferCorr_.max(-0.99);

    //Td *= 1.0/max(1.0,1.001*mag(heatTransferCorr_));
    driftTemp_ = Td;

    Info << "driftTemp:" << nl
         << "    max(driftTemp)      = " << sum(alpha1*Td).value()
            /sum(alpha1).value() << nl
            << "    max(Tslip)      = " << sum(alpha1*TSlip).value()
               /sum(alpha1).value() << nl
         << "    mean(heatCorr)      = " << sum(alpha1*heatTransferCorr_).value()
            /sum(alpha1).value()
         << "    mean(heatCorrecht)      = " << sum(-alpha1*Td/TSlipM).value()
            /sum(alpha1).value()
         << endl;

    // multiply drift Temperature by heatTransfer coefficient
    return
        //pos(heatTransferCorr_)*
        heatTransfer.K()*
        Td;
}

bool Foam::driftTemperatureModel::writeData(Ostream& os) const
{
    return os.good();
}

// ************************************************************************* //

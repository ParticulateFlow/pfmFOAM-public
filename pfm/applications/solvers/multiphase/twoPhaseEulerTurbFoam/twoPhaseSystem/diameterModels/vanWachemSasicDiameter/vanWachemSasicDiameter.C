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

\*---------------------------------------------------------------------------*/

#include "vanWachemSasicDiameter.H"
#include "addToRunTimeSelectionTable.H"
#include "twoPhaseSystem.H"
#include "PhaseCompressibleTurbulenceModel.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
    defineTypeNameAndDebug(vanWachemSasic, 0);

    addToRunTimeSelectionTable
    (
        diameterModel,
        vanWachemSasic,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::vanWachemSasic::vanWachemSasic
(
    const dictionary& diameterProperties,
    const phaseModel& phase
)
:
    diameterModel(diameterProperties, phase),
    d0_("d", dimLength, diameterProperties_),
    H_("H", dimensionSet(1, 2, -2, 0, 0), diameterProperties_),
    delta_("delta", dimLength, diameterProperties_),
    rhoB_("rhoB", dimensionSet(1, -3, 0, 0, 0), diameterProperties_),
    g_("g", dimensionSet(0, 1, -2, 0, 0), diameterProperties_),
    nu_("nu", dimensionSet(0, 0, 0, 0, 0), diameterProperties_),
    E_("E", dimensionSet(1, -1, -2, 0, 0), diameterProperties_),
    d_
    (
        IOobject
        (
            IOobject::groupName("d", phase.name()),
            phase_.U().time().timeName(),
            phase_.U().mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phase_.U().mesh(),
        d0_
     )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diameterModels::vanWachemSasic::~vanWachemSasic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::diameterModels::vanWachemSasic::correct()
{
    const twoPhaseSystem& fluid = refCast<const twoPhaseSystem>(phase_.fluid());
    const volScalarField& rho2 = fluid.otherPhase(phase_).rho();
    // cont. Phase velocity
    const volVectorField& Uc = fluid.otherPhase(phase_).U();
    // const volScalarField& beta = phase_.U().mesh().lookupObject<volScalarField>("Kd");
    volScalarField alpha1(phase_);
    alpha1.max(0);
    volScalarField alpha2(1.0-alpha1);
    
    volScalarField cellVolume
    (
        IOobject
        (
            "cellVolume",
            phase_.U().mesh().time().timeName(),
            phase_.U().mesh()
        ),
        phase_.U().mesh(),
        dimensionedScalar("one", dimLength*dimLength*dimLength, 1)
    );
    cellVolume.ref() = phase_.U().mesh().V();
    
    // slip velocity
    volScalarField uSlip(mag(Uc - phase_.U()));
    
    volScalarField a1(0.5236*(rhoB_-rho2)*g_);
    dimensionedScalar k((1.0 - sqr(nu_))/(3.14*E_));
    volScalarField Us(sqrt(phase_.turbulence().k()));
    volScalarField a2a(-0.166*pow(3.14*pow(Us,6)*pow(rhoB_,3)/sqr(k),0.2));
    volScalarField a2b(-0.1728*rho2*sqr(uSlip)*pow(alpha2,-4.8));
    
    // Apply drag correction
    if (phase_.U().mesh().foundObject<volScalarField>("dragCorr")) {
        const volScalarField& dragCorr = phase_.U().mesh().lookupObject<volScalarField>("dragCorr");
        a2b *= (scalar(1.0) - dragCorr);
    }
    volScalarField a2(a2a+a2b);
    
    dimensionedScalar a3(H_/(24.*sqr(delta_)));
    
    volScalarField det(sqr(a2) - 4.0*a1*a3);
    volScalarField posDet(pos(det));
    det.max(0);
    
    volScalarField d
    (
        posDet
       *(
            - a2
            - sqrt(det)
        )
       /(2.0*a1)
    );
    d = min(0.33*cbrt(cellVolume),d);
    d.max(d0_.value());
    
    // underrelaxation of diameter
    d_ = 0.9*d_ + 0.1*d;
    d_.correctBoundaryConditions();
    Info << "vanWachem & Sasic diameter model: max(d) = " << max(d_).value() << ", min(d) = " << min(d_).value() << endl;
}

bool Foam::diameterModels::vanWachemSasic::read(const dictionary& phaseProperties)
{
    diameterModel::read(phaseProperties);

    diameterProperties_.lookup("d") >> d0_;
    diameterProperties_.lookup("H") >> H_;
    diameterProperties_.lookup("delta") >> delta_;
    diameterProperties_.lookup("rhoB") >> rhoB_;
    diameterProperties_.lookup("g") >> g_;
    diameterProperties_.lookup("nu") >> nu_;
    diameterProperties_.lookup("E") >> E_;

    return true;
}


// ************************************************************************* //

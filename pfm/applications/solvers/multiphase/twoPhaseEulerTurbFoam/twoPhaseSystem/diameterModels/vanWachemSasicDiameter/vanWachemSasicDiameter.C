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
/*
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
    
    d_ = 1.01*d0_;
    
    const cellList& cells = phase_.U().mesh().cells();
    forAll(cells,cellI)
    {
        
        // Newton Method to find root in the interval [d_particle,cellSize]
        // initial conditions
        scalar dx = 100.;
        int iter = 0;
        // Newton iterations
        while ((mag(dx) > 5.0e-6) && (iter < 20) && ((d0_.value() - d_[cellI])*(d_[cellI] - cbrt(cellVolume[cellI]))>0.0)) {
            scalar f  = a1[cellI]*sqr(d_[cellI]) + a2[cellI]*d_[cellI] + a3.value();
            scalar df = 2.0*a1[cellI]*d_[cellI] + a2[cellI];
            dx = -f/df;
            d_[cellI] += dx;
            iter++;
            // Info << "vanWachem & Sasic diameter model: dx = " << dx << ", iter:" << iter  << ", d: " << d_[cellI] << ", alpha: " << alpha1[cellI] << endl;
            // Info << "rho2 = " << rho2[cellI] << ", Uc:" << mag(Uc[cellI])  << ", Us: " << mag(Us[cellI]) << endl;
            // Info << "a1 = " << a1[cellI]*pow(d_[cellI],4) << ", a2:" << a2[cellI]*pow(d_[cellI],3)  << ", a3: " << (a3.value()*pow(d_[cellI],2)) << ", a4: " << a4[cellI] << endl;
            // Info << "a1 = " << a1[cellI] << ", a2:" << a2[cellI]  << ", a3: " << (a3.value()) << ", k: " << k.value() << ", E: " << E_.value() << ", nu: " << nu_.value() << ", f: " << f << ", df: " << df << endl;
        }
        if ((d0_.value() - d_[cellI])*(d_[cellI] - cbrt(cellVolume[cellI]))<0.0) {
            d_[cellI] = d0_.value();
        }
    }
    // d_.max(d0_.value());
    // d_.min(10*d0_.value());
    Info << "vanWachem & Sasic diameter model: max(d) = " << max(d_).value() << ", min(d) = " << min(d_).value() << endl;
}
*/
void Foam::diameterModels::vanWachemSasic::correct()
{
    const twoPhaseSystem& fluid = refCast<const twoPhaseSystem>(phase_.fluid());
    const volScalarField& rho2 = fluid.otherPhase(phase_).rho();
    // cont. Phase velocity
    const volVectorField& Uc = fluid.otherPhase(phase_).U().oldTime();
    const volVectorField& Up = phase_.U().oldTime();

    volScalarField alpha1(phase_.oldTime());
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
    volScalarField uSlip(mag(Uc - Up));
    
    // get strain rate
    volSymmTensorField D(dev(symm(fvc::grad(Up))));
    
    volScalarField a1(0.5236*(rhoB_-rho2)*g_);
    dimensionedScalar k((1.0 - sqr(nu_))/(3.14*E_));
    volScalarField a2a(-0.166*pow(3.14*pow((D&&D)*rhoB_,3)/sqr(k),0.2));
    volScalarField a2b(-0.1728*rho2*sqr(uSlip)*pow(alpha2,-4.8));
    
    // Apply drag correction
    if (phase_.U().mesh().foundObject<volScalarField>("dragCorr")) {
        const volScalarField& dragCorr = phase_.U().mesh().lookupObject<volScalarField>("dragCorr");
        a2b *= (scalar(1.0) - dragCorr);
    }
    
    dimensionedScalar a3(H_/(24.*sqr(delta_)));
    
    d_ = 1.01*d0_;
    
    const cellList& cells = phase_.U().mesh().cells();
    forAll(cells,cellI)
    {
        
        // Newton Method to find root in the interval [d_particle,cellSize]
        // initial conditions
        scalar dx = 100.;
        int iter = 0;
        // Newton iterations
        while ((mag(dx) > 5.0e-6) && (iter < 20) && ((d0_.value() - d_[cellI])*(d_[cellI] - 0.33*cbrt(cellVolume[cellI]))>0.0)) {
            scalar f  = a1[cellI]*sqr(d_[cellI]) + a2a[cellI]*pow(d_[cellI],2.2) + a2b[cellI]*d_[cellI] + a3.value();
            scalar df = 2.0*a1[cellI]*d_[cellI] + 2.2*a2a[cellI]*pow(d_[cellI],1.2) + a2b[cellI];
            dx = -f/df;
            d_[cellI] += dx;
            iter++;
            // Info << "vanWachem & Sasic diameter model: dx = " << dx << ", iter:" << iter  << ", d: " << d_[cellI] << ", alpha: " << alpha1[cellI] << endl;
            // Info << "rho2 = " << rho2[cellI] << ", Uc:" << mag(Uc[cellI])  << ", Us: " << mag(Us[cellI]) << endl;
            // Info << "a1 = " << a1[cellI]*pow(d_[cellI],4) << ", a2:" << a2[cellI]*pow(d_[cellI],3)  << ", a3: " << (a3.value()*pow(d_[cellI],2)) << ", a4: " << a4[cellI] << endl;
            // Info << "a1 = " << a1[cellI] << ", a2:" << a2[cellI]  << ", a3: " << (a3.value()) << ", k: " << k.value() << ", E: " << E_.value() << ", nu: " << nu_.value() << ", f: " << f << ", df: " << df << endl;
        }
        if ((d0_.value() - d_[cellI])*(d_[cellI] - 0.33*cbrt(cellVolume[cellI]))<0.0) {
            d_[cellI] = d0_.value();
        }
    }
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

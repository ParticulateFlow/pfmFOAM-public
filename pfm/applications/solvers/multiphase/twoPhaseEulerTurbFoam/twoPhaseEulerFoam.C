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

Application
    twoPhaseEulerFoam

Description
    Solver for a system of 2 compressible fluid phases with one phase
    dispersed, e.g. gas bubbles in a liquid including heat-transfer.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "twoPhaseSystem.H"
#include "PhaseCompressibleTurbulenceModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "fixedValueFvsPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "createTimeControls.H"
    #include "CourantNos.H"
    #include "setInitialDeltaT.H"

    Switch faceMomentum
    (
        pimple.dict().lookupOrDefault<Switch>("faceMomentum", false)
    );

    Switch implicitPhasePressure
    (
        mesh.solverDict(alpha1.name()).lookupOrDefault<Switch>
        (
            "implicitPhasePressure", false
        )
    );
    
    // switch for periodic box simulations
    // handling of p_rgh in case of periodic boxes is wrong,
    // thus, we introduced a second gravity term in the momentum equations
    // in case of periodic box simulations the gravity in the
    // "constant/g" dictionary has to be set to 0, to remove it from the
    // pressure equation!
    Switch periodicBox
    (
        pimple.dict().lookupOrDefault<Switch>("periodicBox", false)
    );
    
    if (!periodicBox) gN *= 0.;
    
    Switch energyEqn
    (
        pimple.dict().lookupOrDefault<Switch>("energyEqn", false)
    );

    #include "pUf/createDDtU.H"
    #include "pU/createDDtU.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNos.H"
        #include "setDeltaT.H"

        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            fluid.solve();
            fluid.correct();

            #include "contErrs.H"

            if (faceMomentum)
            {
                #include "pUf/UEqns.H"
                if (energyEqn) {
                    #include "EEqns.H"
                }
                #include "pUf/pEqn.H"
                #include "pUf/DDtU.H"
            }
            else
            {
                #include "pU/UEqns.H"
                if (energyEqn) {
                    #include "EEqns.H"
                }
                #include "pU/pEqn.H"
                #include "pU/DDtU.H"
            }

            if (pimple.turbCorr())
            {
                fluid.correctTurbulence();
            }
        }
        
        ///////////---------------POST_PROCESS-----------//////////////////////////
        Info<< "particle_ENSTROPHY: "
            << fvc::domainIntegrate( 0.5*magSqr(fvc::curl(U1))).value()
            << endl;
        
        Info<< "air_ENSTROPHY: "
            << fvc::domainIntegrate(0.5*magSqr(fvc::curl(U2))).value()
            << endl;
        if (periodicBox) {
            Info<< "slip_velocity: "
                << - ((
                        fvc::domainIntegrate(alpha2*(U2&gN)).value()
                       /fvc::domainIntegrate(alpha2*mag(gN)).value()
                     )
                   - (
                        fvc::domainIntegrate(alpha1*(U1&gN)).value()
                       /fvc::domainIntegrate(alpha1*mag(gN)).value()
                     ))
                << endl;
        }
        
        Info<< "total momentum: "
            << mag(fvc::domainIntegrate(alpha1*rho1*U1 + alpha2*rho2*U2).value())
            << endl;
        
        Info<< "total solids mass: "
            << mag(fvc::domainIntegrate(alpha1*rho1).value())
            << endl;

        #include "write.H"

        Info<< "ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n\n" << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

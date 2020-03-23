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
    interFoam

Description
    Solver for 2 incompressible, isothermal immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing approach,
    with optional mesh motion and mesh topology changes including adaptive
    re-meshing.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"
#include "wallDist.H"
#include "simpleTestFilter.H"
#include "simpleFilter.H"
#include "laplaceFilter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createAlphaFluxes.H"
    #include "initCorrectPhi.H"
    #include "createUfIfPresent.H"

    turbulence->validate();

    if (!LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "CourantNo.H"
            #include "alphaCourantNo.H"
            #include "setDeltaT.H"
        }

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                mesh.update();

                if (mesh.changing())
                {
                    // Do not apply previous time-step mesh compression flux
                    // if the mesh topology changed
                    if (mesh.topoChanging())
                    {
                        talphaPhi1Corr0.clear();
                    }

                    gh = (g & mesh.C()) - ghRef;
                    ghf = (g & mesh.Cf()) - ghRef;

                    MRF.update();

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        #include "correctPhi.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);

                        mixture.correct();
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            #include "alphaControls.H"
            #include "alphaEqnSubCycle.H"

            mixture.correct();

            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                volScalarField nu(mesh.lookupObject<volScalarField>("nu"));
                
                tmp<volTensorField> tgradU(fvc::grad(U));
                const volTensorField& gradU(tgradU());
                volSymmTensorField D(symm(gradU));
                 // Dynamic adjustment of Cmu
                /*
                volScalarField Lij  = filter_(magSqr(U)) - magSqr(filter_(U));
                Lij.max(0);
                volScalarField Mij = sqr(deltaF_)*(4.0*magSqr(filter_(D)) - filter_(magSqr(D)));
                volScalarField MijMij = filterS_(sqr(Mij));
                MijMij.max(SMALL);
        
                volScalarField CmuT = 0.5*(filterS_(Lij * Mij)/(MijMij));
                CmuT.min(1.0);
                CmuT.max(0.001);
                
                Cmu_ = sqrt(CmuT);
            */
                volScalarField nutSigmaCorr = -sqr(Csigma_)*(mixture.sigmaK())*(fvc::laplacian(alpha1))/rho;
                nutSigmaCorr.max(SMALL);
                
                volScalarField a(Ceps_/deltaF_);
                volScalarField b((2.0/3.0)*tr(D));
                volScalarField c(Cmu_*deltaF_*(2*(dev(D) && D) + nutSigmaCorr*mixture.nearInterface()));
                
                volScalarField k(sqr((-b + sqrt(sqr(b) + 4*a*c))/(2*a)));
                
                // standard
                nutSigma_ =  Cmu_*deltaF_*sqrt(k);
                                 
                
                // WALE
                /*
                volSymmTensorField Sijd(dev(symm(gradU&gradU)));
                volScalarField strain(pow3(magSqr(Sijd))
                                    /sqr(
                                        pow(magSqr(dev(D)),5.0/2.0)
                                      + pow(magSqr(Sijd),5.0/4.0)
                                      + dimensionedScalar
                                           (
                                               "small",
                                               dimensionSet(0, 0, -5, 0, 0),
                                               small
                                           )
                                      )
                                    );
                nutSigma_ =  sqr(Cmu_*deltaF_)
                           * sqrt(
                                     2.0*strain
                                   + nutSigmaCorr
                                );
                */
                nutSigma_.correctBoundaryConditions();
                // Limit nut
                nutSigma_ = min(nutSigma_,1.0e5*nu);
                nutSigma_.max(SMALL); 
 
                // Dynamic adjustment of Cst
                
                volSymmTensorField Dhat(filter_(D));
                volScalarField sigmaKhat(filter_(mixture.sigmaK()));
                volScalarField nutSigmaCorrH = -sqr(Csigma_)*sigmaKhat*filter_(fvc::laplacian(alpha1)/rho);
                nutSigmaCorrH.max(SMALL);
                
                volScalarField aH(Ceps_/(2.0*deltaF_));
                volScalarField bH((2.0/3.0)*tr(Dhat));
                volScalarField cH(2.0*Cmu_*deltaF_*(2.0*(dev(Dhat) && Dhat) + nutSigmaCorrH*mixture.nearInterface()));
                
                volScalarField kH(sqr((-bH + sqrt(sqr(bH) + 4*aH*cH))/(2*aH)));
                
                // standard
                volScalarField nutSigmaHat =  2.0*Cmu_*deltaF_*sqrt(kH);
                
               
                nutSigmaHat = min(nutSigmaHat,1.0e5*nu);
                nutSigmaHat.max(SMALL);
                
                volVectorField gradAlpha = fvc::grad(alpha1);
                volVectorField gradAlphaHat = filter_(gradAlpha);
 
                volVectorField MijS = sigmaKhat*gradAlphaHat*sqrt(nutSigmaHat/nu)
                                    - filter_(mixture.sigmaK()*gradAlpha*sqrt(nutSigma_/nu));
 
                volVectorField LijS = 2.0*(filter_(mixture.sigmaK()*gradAlpha) - sigmaKhat*gradAlphaHat);
                volScalarField MijMijS = filterS_(MijS&MijS);
                MijMijS.max(ROOTVSMALL);
              
                Cst_ = mixture.nearInterface()*(filterS_(LijS&MijS))/MijMijS;
                Cst_ = 0.5*(Cst_+mag(Cst_));
            
                Info << "max(nut) = " << max(nutSigma_).value() << nl
                     << "min(nut) = " << min(nutSigma_).value() << endl;
                
                corrSurfaceTensionForce_ = (
                                               scalar(1.0)
                                             + Cst_*sqrt(nutSigma_/nu)
                                              *mixture.nearInterface()
                                           );
                corrSurfaceTensionForce_.max(0.01);
                corrSurfaceTensionForce_.min(100.0);
                turbulence->correct();
            }
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

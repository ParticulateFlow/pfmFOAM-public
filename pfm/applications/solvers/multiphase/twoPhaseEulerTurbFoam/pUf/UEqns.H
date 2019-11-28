Info<< "Constructing face momentum equations" << endl;

MRF.correctBoundaryVelocity(U1);
MRF.correctBoundaryVelocity(U2);
MRF.correctBoundaryVelocity(U);

fvVectorMatrix U1Eqn(U1, rho1.dimensions()*U1.dimensions()*dimVol/dimTime);
fvVectorMatrix U2Eqn(U2, rho2.dimensions()*U2.dimensions()*dimVol/dimTime);

// modification for periodic box simulations
volScalarField rhoMix = 0.0*rho1;
volVectorField totalMomentum = 0.0*U1;
dimensionedScalar deltaT = U1.mesh().time().deltaT();
if (periodicBox) {
    rhoMix = fvc::domainIntegrate(alpha1*rho1+alpha2*rho2)
                / fvc::domainIntegrate(unity);
    
    totalMomentum = fvc::domainIntegrate(alpha1*rho1*U1 + alpha2*rho2*U2)
                /fvc::domainIntegrate(alpha1*rho1+alpha2*rho2);
}

volVectorField uSlipV(U2 - U1);

{
    volScalarField Vm(fluid.Vm());
    
    fvVectorMatrix UgradU1
    (
        fvm::div(phi1, U1) - fvm::Sp(fvc::div(phi1), U1)
      + MRF.DDt(U1)
    );

    fvVectorMatrix UgradU2
    (
        fvm::div(phi2, U2) - fvm::Sp(fvc::div(phi2), U2)
      + MRF.DDt(U2)
    );

    {
        U1Eqn =
        (
            fvm::div(alphaRhoPhi1, U1)
          + fvm::SuSp(-fvc::div(alphaRhoPhi1), U1)
          + MRF.DDt(alpha1*rho1, U1)
          + phase1.turbulence().divDevRhoReff(U1)
          + Vm*(UgradU1 - (UgradU2 & U2))
          // own gravity for periodic box
          // (gN = 0, if pimple switch "periodicBox" is false)
          - alpha1*(rho1-rhoMix)*gN
          // in case of periodic box, adjust CG-velocity
          + alpha1*rho1*totalMomentum/deltaT
          - fvOptions(alpha1, rho1, U1)
        );
        U1Eqn.relax();
        fvOptions.constrain(U1Eqn);
        U1.correctBoundaryConditions();
        fvOptions.correct(U1);
    }

    {
        U2Eqn =
        (
            fvm::div(alphaRhoPhi2, U2)
          + fvm::SuSp(-fvc::div(alphaRhoPhi2), U2)
          + MRF.DDt(alpha2*rho2, U2)
          + phase2.turbulence().divDevRhoReff(U2)
          + Vm*(UgradU2 - (UgradU1 & U1))     
          // own gravity for periodic box
          // (gN = 0, if pimple switch "periodicBox" is false)
          - alpha2*(rho2-rhoMix)*gN
          // in case of periodic box, adjust CG-velocity
          + alpha2*rho2*totalMomentum/deltaT
          - fvOptions(alpha2, rho2, U2)
        );
        U2Eqn.relax();
        fvOptions.constrain(U2Eqn);
        U2.correctBoundaryConditions();
        fvOptions.correct(U2);
    }
}

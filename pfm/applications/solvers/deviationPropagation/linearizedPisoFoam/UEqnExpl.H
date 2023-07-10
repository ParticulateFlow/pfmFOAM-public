// Solve the Momentum equation

fvVectorMatrix UEqn
(
    fvm::ddt(U) + fvc::div(phi_ts, U)
    + fvc::div(phi,U_ts)
    - fvc::laplacian(turbulence->nu(), U)
 ==
    fvOptions(U)
);

UEqn.relax();

fvOptions.constrain(UEqn);

volTensorField gKuu = fvc::grad(K_uu);
volVectorField KgradU = K_uu & fvc::grad(U_ts);

fvVectorMatrix KuuEqn
(
    fvm::ddt(K_uu)
    + fvc::div(phi_ts, K_uu)
    + KgradU
    - fvc::laplacian(turbulence->nu(), K_uu)
 ==
    fvOptions(K_uu)
);

KuuEqn.relax();

solve(KuuEqn);

fvOptions.correct(K_uu);

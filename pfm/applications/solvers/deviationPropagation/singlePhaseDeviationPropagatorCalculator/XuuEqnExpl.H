volTensorField gXuu = fvc::grad(X_uu);
volVectorField XgradU = X_uu & fvc::grad(U_ts);

fvVectorMatrix XuuEqn
(
    fvm::ddt(X_uu)
    + fvc::div(phi_ts, X_uu)
    + XgradU
    - fvc::laplacian(turbulence->nu(), X_uu)
 ==
    fvOptions(X_uu)
);

XuuEqn.relax();

solve(XuuEqn);
//solve(XuuEqn == -fvc::grad(X_pu));


fvOptions.correct(X_uu);

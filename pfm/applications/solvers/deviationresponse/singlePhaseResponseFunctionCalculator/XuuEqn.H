fvVectorMatrix XuuEqn
(
    fvm::ddt(X_uu)
    + fvm::div(phi_ts, X_uu)
 //   + turbulence->divDevReff(X_uu)
    - fvm::laplacian(turbulence->nu(), X_uu)
 ==
    fvOptions(X_uu)
);

XuuEqn.relax();

tmp<volTensorField> tDiagXuu = tensor(I)*XuuEqn.A();
tDiagXuu.ref() += fvc::grad(U_ts);
tmp<volTensorField> trDiagXuu = inv(tDiagXuu());
trDiagXuu.ref().rename("rAU");

fvOptions.constrain(XuuEqn);

volVectorField gradXpu(fvc::grad(X_pu));

for (int XuuCorr=0; XuuCorr<nXuuCorr; XuuCorr++)
{
    X_uu = trDiagXuu() & (XuuEqn.H() - gradXpu);
}
X_uu.correctBoundaryConditions();

fvOptions.correct(X_uu);

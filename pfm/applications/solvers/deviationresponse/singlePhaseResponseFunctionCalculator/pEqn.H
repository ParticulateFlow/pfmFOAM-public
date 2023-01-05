volTensorField gXuu = fvc::grad(X_uu);
volTensorField gUdb = fvc::grad(U_db);

RHS = gXuu.component(tensor::XX)*gUdb.component(tensor::XX)
    + gXuu.component(tensor::XY)*gUdb.component(tensor::YX)
    + gXuu.component(tensor::XZ)*gUdb.component(tensor::ZX)
    + gXuu.component(tensor::YX)*gUdb.component(tensor::XY)
    + gXuu.component(tensor::YY)*gUdb.component(tensor::YY)
    + gXuu.component(tensor::YZ)*gUdb.component(tensor::ZY)
    + gXuu.component(tensor::ZX)*gUdb.component(tensor::XZ)
    + gXuu.component(tensor::ZY)*gUdb.component(tensor::YZ)
    + gXuu.component(tensor::ZZ)*gUdb.component(tensor::ZZ)
    + gUdb.component(tensor::XX)*gXuu.component(tensor::XX)
    + gUdb.component(tensor::XY)*gXuu.component(tensor::YX)
    + gUdb.component(tensor::XZ)*gXuu.component(tensor::ZX)
    + gUdb.component(tensor::YX)*gXuu.component(tensor::XY)
    + gUdb.component(tensor::YY)*gXuu.component(tensor::YY)
    + gUdb.component(tensor::YZ)*gXuu.component(tensor::ZY)
    + gUdb.component(tensor::ZX)*gXuu.component(tensor::XZ)
    + gUdb.component(tensor::ZY)*gXuu.component(tensor::YZ)
    + gUdb.component(tensor::ZZ)*gXuu.component(tensor::ZZ);


tmp<fvScalarMatrix> tXpuEqn
(
    fvm::laplacian(X_pu) == - RHS
);
fvScalarMatrix XpuEqn = tXpuEqn.ref();

XpuEqn.setReference(pRefCell, pRefValue);
XpuEqn.solve();

X_uu = trDiagXuu() & (XuuEqn.H() - fvc::grad(X_pu));
X_uu.correctBoundaryConditions();

Info<< "max(X_uu) = " << max(X_uu).value() << endl;

    
/*
// standard PISO pressure treatment
tmp<volVectorField> tHbyA = constrainHbyA(trDiagXuu()& (XuuEqn.H() - X_uu.oldTime() / deltaT), X_uu, X_pu);
volVectorField& HbyA = tHbyA.ref();
surfaceScalarField phiHbyA("phiHbyA",fvc::flux(HbyA));

// Non-orthogonal pressure corrector loop
while (piso.correctNonOrthogonal())
{
    // Pressure corrector

    fvScalarMatrix pEqn
    (
        fvm::laplacian(trDiagXuu(), X_pu) == fvc::div(phiHbyA)
    );

    pEqn.setReference(pRefCell, pRefValue);

    pEqn.solve(mesh.solver(X_pu.select(piso.finalInnerIter())));

 //   if (piso.finalNonOrthogonalIter())
 //   {
 //       phiX = phiHbyA - pEqn.flux();
 //   }
}

X_uu = HbyA - (trDiagXuu()&fvc::grad(X_pu)) + (trDiagXuu()&X_uu.oldTime() / deltaT);
X_uu.correctBoundaryConditions();
fvOptions.correct(X_uu);
*/

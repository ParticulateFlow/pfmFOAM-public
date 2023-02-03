// standard PISO pressure treatment
tmp<volVectorField> tHbyA = constrainHbyA(trDiagXuu()& XuuEqn.H(), X_uu, X_pu);
volVectorField& HbyA = tHbyA.ref();
surfaceScalarField phiHbyA("phiHbyA",fvc::flux(HbyA));

// Non-orthogonal pressure corrector loop
while (piso.correctNonOrthogonal())
{
    // Pressure corrector

    fvScalarMatrix XpuEqn
    (
        fvm::laplacian(trDiagXuu(), X_pu) == fvc::div(phiHbyA)
    );

    XpuEqn.setReference(pRefCell, pRefValue);

    XpuEqn.solve(mesh.solver(X_pu.select(piso.finalInnerIter())));

 //   if (piso.finalNonOrthogonalIter())
 //   {
 //       phiX = phiHbyA - pEqn.flux();
 //   }
}

X_uu = HbyA - (trDiagXuu()&fvc::grad(X_pu));
X_uu.correctBoundaryConditions();
fvOptions.correct(X_uu);

Info<< "max(X_uu) = " << max(X_uu).value() << endl;

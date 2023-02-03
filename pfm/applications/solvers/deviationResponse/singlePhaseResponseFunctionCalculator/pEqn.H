/*
// different discretization as in momentum equation; would it be better to have div(div(...))
volTensorField gXuu = fvc::grad(X_uu);
volTensorField gUts = fvc::grad(U_ts);

RHS = gXuu.component(tensor::XX)*gUts.component(tensor::XX)
    + gXuu.component(tensor::XY)*gUts.component(tensor::YX)
    + gXuu.component(tensor::XZ)*gUts.component(tensor::ZX)
    + gXuu.component(tensor::YX)*gUts.component(tensor::XY)
    + gXuu.component(tensor::YY)*gUts.component(tensor::YY)
    + gXuu.component(tensor::YZ)*gUts.component(tensor::ZY)
    + gXuu.component(tensor::ZX)*gUts.component(tensor::XZ)
    + gXuu.component(tensor::ZY)*gUts.component(tensor::YZ)
    + gXuu.component(tensor::ZZ)*gUts.component(tensor::ZZ)
    + gUts.component(tensor::XX)*gXuu.component(tensor::XX)
    + gUts.component(tensor::XY)*gXuu.component(tensor::YX)
    + gUts.component(tensor::XZ)*gXuu.component(tensor::ZX)
    + gUts.component(tensor::YX)*gXuu.component(tensor::XY)
    + gUts.component(tensor::YY)*gXuu.component(tensor::YY)
    + gUts.component(tensor::YZ)*gXuu.component(tensor::ZY)
    + gUts.component(tensor::ZX)*gXuu.component(tensor::XZ)
    + gUts.component(tensor::ZY)*gXuu.component(tensor::YZ)
    + gUts.component(tensor::ZZ)*gXuu.component(tensor::ZZ);

dimensionedScalar rdt("zero",dimensionSet(0,0,-1,0,0),1.0/deltaT.value());
// is the second term on the right-hand side necessary? shouldn't make a difference for response to initial state (div free per definition), but for response to time-derivative
tmp<fvScalarMatrix> tXpuEqn
(
    fvm::laplacian(X_pu) == - RHS //+ rdt*fvc::div(X_uu)
);
fvScalarMatrix XpuEqn = tXpuEqn.ref();

XpuEqn.setReference(pRefCell, pRefValue);
XpuEqn.solve();

//X_uu = trDiagXuu() & (XuuEqn.H() - fvc::grad(X_pu));


gradXpu = fvc::grad(X_pu);
for (int XuuCorr=0; XuuCorr<nXuuCorr; XuuCorr++)
{
    X_uu = trDiagXuu() & (XuuEqn.H() - gradXpu);
    X_uu.correctBoundaryConditions();
}


X_uu.correctBoundaryConditions();

Info<< "max(X_uu) = " << max(X_uu).value() << endl;

*/
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
Info<< "max(X_uu) = " << max(X_uu).value() << endl;

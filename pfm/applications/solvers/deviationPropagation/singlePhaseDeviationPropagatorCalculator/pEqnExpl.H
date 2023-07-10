// different discretization as in momentum equation; would it be better to have div(div(...))

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

// Non-orthogonal pressure corrector loop
while (piso.correctNonOrthogonal())
{
    tmp<fvScalarMatrix> tXpuEqn
    (
        fvm::laplacian(X_pu) == - RHS
    );
    fvScalarMatrix XpuEqn = tXpuEqn.ref();

    XpuEqn.setReference(pRefCell, pRefValue);
    XpuEqn.solve();
}

X_uu -= 1.0/rdt*fvc::grad(X_pu);
X_uu.correctBoundaryConditions();

Info<< "max(X_uu) = " << max(X_uu).value() << endl;

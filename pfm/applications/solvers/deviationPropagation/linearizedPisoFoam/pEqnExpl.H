// different discretization as in momentum equation; would it be better to have div(div(...))

volTensorField gUts = fvc::grad(U_ts);

RHS = gKuu.component(tensor::XX)*gUts.component(tensor::XX)
    + gKuu.component(tensor::XY)*gUts.component(tensor::YX)
    + gKuu.component(tensor::XZ)*gUts.component(tensor::ZX)
    + gKuu.component(tensor::YX)*gUts.component(tensor::XY)
    + gKuu.component(tensor::YY)*gUts.component(tensor::YY)
    + gKuu.component(tensor::YZ)*gUts.component(tensor::ZY)
    + gKuu.component(tensor::ZX)*gUts.component(tensor::XZ)
    + gKuu.component(tensor::ZY)*gUts.component(tensor::YZ)
    + gKuu.component(tensor::ZZ)*gUts.component(tensor::ZZ)
    + gUts.component(tensor::XX)*gKuu.component(tensor::XX)
    + gUts.component(tensor::XY)*gKuu.component(tensor::YX)
    + gUts.component(tensor::XZ)*gKuu.component(tensor::ZX)
    + gUts.component(tensor::YX)*gKuu.component(tensor::XY)
    + gUts.component(tensor::YY)*gKuu.component(tensor::YY)
    + gUts.component(tensor::YZ)*gKuu.component(tensor::ZY)
    + gUts.component(tensor::ZX)*gKuu.component(tensor::XZ)
    + gUts.component(tensor::ZY)*gKuu.component(tensor::YZ)
    + gUts.component(tensor::ZZ)*gKuu.component(tensor::ZZ);

dimensionedScalar rdt("zero",dimensionSet(0,0,-1,0,0),1.0/deltaT.value());
tmp<fvScalarMatrix> tKpuEqn
(
    fvm::laplacian(K_pu) == - RHS
);
fvScalarMatrix KpuEqn = tKpuEqn.ref();

KpuEqn.setReference(pRefCell, pRefValue);
KpuEqn.solve();

K_uu -= 1.0/rdt*fvc::grad(K_pu);
K_uu.correctBoundaryConditions();

Info<< "max(K_uu) = " << max(K_uu).value() << endl;

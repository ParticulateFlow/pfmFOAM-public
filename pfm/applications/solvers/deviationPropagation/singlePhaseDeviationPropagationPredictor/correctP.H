dimensionedScalar dummyT("one",dimensionSet(0,0,1,0,0,0,0),1.0);

fvScalarMatrix pEqn
(
    fvm::laplacian(p) == fvc::div(U)/dummyT
);

pEqn.setReference(pRefCell, pRefValue);
pEqn.solve();

U -= dummyT*fvc::grad(p);
U.correctBoundaryConditions();

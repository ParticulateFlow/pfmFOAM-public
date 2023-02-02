volTensorField X_uu_integrated
(
    IOobject
    (
        "X_uu_integrated",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedTensor("zero",dimensionSet(0,0,0,0,0,0,0),tensor::zero)
);

volTensorField X_uu_integrated_target
(
    IOobject
    (
        "X_uu_integrated_target",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedTensor("zero",dimensionSet(0,0,0,0,0,0,0),tensor::zero)
);

volTensorField deltaX_uu
(
    IOobject
    (
        "deltaX_uu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedTensor("zero",dimensionSet(0,0,-3,0,0,0,0),tensor::zero)
);

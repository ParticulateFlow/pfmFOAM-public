volTensorField K_uu_integrated
(
    IOobject
    (
        "K_uu_integrated",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedTensor("zero",dimensionSet(0,0,0,0,0,0,0),tensor::zero)
);

volTensorField K_uu_integrated_target
(
    IOobject
    (
        "K_uu_integrated_target",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedTensor("zero",dimensionSet(0,0,0,0,0,0,0),tensor::zero)
);

IOdictionary dbProperties
(
    IOobject
    (
        "dataBaseProperties",
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

scalar minKuu(dbProperties.lookupOrDefault<scalar>("minKuu",1e-5));

    IOdictionary fieldDifferenceProperties
    (
        IOobject
        (
            "fieldDifferenceProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    wordList fieldNames(fieldDifferenceProperties.lookupOrDefault<wordList>("fieldNames",wordList(0)));
    word fieldType(fieldDifferenceProperties.lookup("fieldType"));

    if (fieldNames.size() != 2)
    {
        FatalError << "Please provide exactly two field names\n" << abort(FatalError);
    }

    if (fieldType != "volScalarField" && fieldType != "volVectorField" && fieldType != "volSymmTensorField")
    {
        FatalError << "Field type has to be either 'volScalerField' or 'volVectorField' or 'volSymmTensorField'\n" << abort(FatalError);
    }

    volScalarField diff 
    (
        IOobject
        (
            "diff",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero",dimensionSet(0,0,0,0,0,0,0), 0.0)
    );

    volScalarField norm0 
    (
        IOobject
        (
            "norm1",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero",dimensionSet(0,0,0,0,0,0,0), 0.0)
    );

    volScalarField norm1
    (
        IOobject
        (
            "norm2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero",dimensionSet(0,0,0,0,0,0,0), 0.0)
    );

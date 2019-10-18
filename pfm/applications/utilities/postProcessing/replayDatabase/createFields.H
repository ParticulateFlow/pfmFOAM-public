    IOdictionary replayDatabaseProperties
    (
        IOobject
        (
            "replayDatabaseProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    label samplestep(replayDatabaseProperties.lookupOrDefault<label>("sampleStep",1));

    scalar startTime(replayDatabaseProperties.lookupOrDefault<scalar>("startTime",0.0));

    scalar endTime(replayDatabaseProperties.lookupOrDefault<scalar>("endTime",0.0));

    bool limitTime(replayDatabaseProperties.lookupOrDefault<bool>("limitTime",false));

    word databaseName(replayDatabaseProperties.lookupOrDefault<word>("databaseName","database"));

    wordList distDataScalarFieldNames_(replayDatabaseProperties.lookup("distDataScalarFields"));

    wordList distDataVectorFieldNames_(replayDatabaseProperties.lookup("distDataVectorFields"));

    PtrList<volScalarField> distDataScalarFields_;

    PtrList<volVectorField> distDataVectorFields_;

    Foam::Time dbTime(fileName(databaseName), "", "../system", "../constant", false);
    instantList timeDirs(dbTime.times());
    dbTime.setTime(timeDirs[0],0);
/*
    for (int i=0; i<distDataScalarFieldNames_.size(); i++)
    {
        distDataScalarFields_.append
        (
            new volScalarField
            (
                IOobject
                (
                    distDataScalarFieldNames_[i],
                    dbTime.timePath(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                 ),
                 mesh
            )
        );
    }

    for (int i=0; i<distDataVectorFieldNames_.size(); i++)
    {
        distDataVectorFields_.append
        (
            new volVectorField
            (
                IOobject
                (
                    distDataVectorFieldNames_[i],
                    dbTime.timePath(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                 ),
                 mesh
            )
        );
    }
*/

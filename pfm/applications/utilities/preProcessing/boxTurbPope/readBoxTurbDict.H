    Info<< "Reading boxTurbDict\n" << endl;

    IOdictionary boxTurbDict
    (
        IOobject
        (
            "boxTurbDict",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    scalar Cp(readScalar(boxTurbDict.lookup("Cp")));

    scalar p0(readScalar(boxTurbDict.lookup("p0")));

    scalar beta(readScalar(boxTurbDict.lookup("beta")));
    
    scalar nu(readScalar(boxTurbDict.lookup("nu")));

    scalar cL(readScalar(boxTurbDict.lookup("cL")));

    scalar cEta(readScalar(boxTurbDict.lookup("cEta")));

    scalar e(readScalar(boxTurbDict.lookup("e")));

    scalar epsilon(readScalar(boxTurbDict.lookup("epsilon")));

    Info << "C    = " << Cp << nl
         << "p0   = " << p0<< nl
         << "beta = " << beta << nl
         << "nu   = " << nu << nl
         << "cL   = " << cL << nl
         << "cEta = " << cEta << nl
         << "e    = " << e << nl
         << "eps  = " << epsilon << nl << endl;

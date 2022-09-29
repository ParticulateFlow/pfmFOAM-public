{
    volScalarField& he1 = thermo1.he();
    volScalarField& he2 = thermo2.he();

    volScalarField Cpv1("Cpv1", thermo1.Cpv());
    volScalarField Cpv2("Cpv2", thermo2.Cpv());

    volScalarField Kh(fluid.Kh());

    volScalarField KhTdrift(fluid.KhTdrift());

    if (mesh.foundObject<volScalarField>("heatTransferCorr")) {
        Info<< "Applying heat transfer correction..." << endl;
        const volScalarField& heatTransferCorr = mesh.lookupObject<volScalarField>("heatTransferCorr");
        Kh *= (scalar(1.0) - heatTransferCorr);
    }

    Kd = fluid.Kd();

    volVectorField KdUdrift(fluid.KdUdrift());

    if (mesh.foundObject<volScalarField>("dragCorr")) {
        Info<< "Applying drag correction..." << endl;
        const volScalarField& dragCorr = mesh.lookupObject<volScalarField>("dragCorr");
        Kd *= (scalar(1.0) - dragCorr);
    }

    volVectorField dragEx = Kd*(U2 - U1);

    const dictionary& thermophase1 = mesh.lookupObject<IOdictionary>
    (
       "thermophysicalProperties." +  phase1.name()
    );

    dimensionedScalar Prt1
    (
        "Prt1",
        dimensionSet(0,0,0,0,0),
        thermophase1.subDict("mixture").subDict("transport").lookupOrDefault<scalar>("Prt",1.0)
    );

    volScalarField alphat1 = phase1.turbulence().mut()/Prt1;

    const dictionary& thermophase2 = mesh.lookupObject<IOdictionary>
    (
       "thermophysicalProperties." +  phase2.name()
    );

    dimensionedScalar Prt2
    (
        "Prt2",
        dimensionSet(0,0,0,0,0),
        thermophase2.subDict("mixture").subDict("transport").lookupOrDefault<scalar>("Prt",1.0)
    );

    volScalarField alphat2 = phase2.turbulence().mut()/Prt2;


    const dictionary& turbulencephase1 = mesh.lookupObject<IOdictionary>
    (
       "turbulenceProperties." +  phase1.name()
    );
    const word turbulence1Groupname = turbulencephase1.lookup("simulationType");
    const word turbulence1name = turbulencephase1.subDict("RAS").lookup("RASModel");


    if (turbulence1Groupname == "RAS" && turbulence1name == "kineticTheory")
    {
       alphat1 = mesh.lookupObject<volScalarField>("alphaturb." + phase1.name());

    }

    volSymmTensorField viscousHeat1 = phase1.turbulence().devRhoReff();
    volSymmTensorField viscousHeat2 = phase2.turbulence().devRhoReff();

    fvScalarMatrix E1Eqn
    (
        fvm::ddt(alpha1, rho1, he1) + fvm::div(alphaRhoPhi1, he1)
      - fvm::Sp(contErr1, he1)

      + fvc::ddt(alpha1, rho1, K1) + fvc::div(alphaRhoPhi1, K1)
      - contErr1*K1
      + (
            he1.name() == thermo1.phasePropertyName("e")
          ? fvc::div(fvc::absolute(alphaPhi1, alpha1, U1), p)
          + p*fvc::ddt(alpha1) + alpha1*p*fvc::div(fvc::absolute(phi1,U1))
          : -alpha1*dpdt + alpha1*p*fvc::div(fvc::absolute(phi1,U1))
        )

      - fvm::laplacian
        (
            fvc::interpolate(alpha1)
           *fvc::interpolate(thermo1.alpha()
                             + alphat1
                             ),
            he1
        )
      + fvc::div(viscousHeat1&U1)
      - (U1&dragEx)
   
    );

    E1Eqn.relax();

    E1Eqn -=
    (
        Kh*(thermo2.T() - thermo1.T())
      + Kh*he1/Cpv1
      - fvm::Sp(Kh/Cpv1, he1)
      + alpha1*rho1*(U1&g)
      + fvOptions(alpha1, rho1, he1)
    );

    fvScalarMatrix E2Eqn
    (
        fvm::ddt(alpha2, rho2, he2) + fvm::div(alphaRhoPhi2, he2)
      - fvm::Sp(contErr2, he2)

      + fvc::ddt(alpha2, rho2, K2) + fvc::div(alphaRhoPhi2, K2)
      - contErr2*K2
      + (
            he2.name() == thermo2.phasePropertyName("e")
          ? fvc::div(fvc::absolute(alphaPhi2, alpha2, U2), p)
          + p*fvc::ddt(alpha2) + alpha2*p*fvc::div(fvc::absolute(phi2,U2))
          : -alpha2*dpdt + alpha2*p*fvc::div(fvc::absolute(phi2,U2))
        )

      - fvm::laplacian
        (
            fvc::interpolate(alpha2)
           *fvc::interpolate(thermo2.alpha()
                             + alphat2
                             ),
            he2
        )

       + fvc::div(viscousHeat2&U2)
       + (U2&dragEx)
    );

    E2Eqn.relax();

    E2Eqn -=
    (
        Kh*(thermo1.T() - thermo2.T())
      + Kh*he2/Cpv2
      - fvm::Sp(Kh/Cpv2, he2)
      + alpha2*rho2*(U2&g)
      + fvOptions(alpha2, rho2, he2)
    );

    fvOptions.constrain(E1Eqn);
    E1Eqn.solve();
    fvOptions.correct(he1);

    fvOptions.constrain(E2Eqn);
    E2Eqn.solve();
    fvOptions.correct(he2);

    thermo1.correct();
    Info<< "min " << thermo1.T().name()
        << " " << min(thermo1.T()).value() << endl;

    thermo2.correct();
    Info<< "min " << thermo2.T().name()
        << " " << min(thermo2.T()).value() << endl;

    Info<< "T_diff: "
        << fvc::domainIntegrate((thermo2.T() - thermo1.T())/(dimensionedScalar("500", dimensionSet(0,0,0,1,0), 500.0) - thermo1.T())).value()
            /fvc::domainIntegrate(unity).value() << endl;

}

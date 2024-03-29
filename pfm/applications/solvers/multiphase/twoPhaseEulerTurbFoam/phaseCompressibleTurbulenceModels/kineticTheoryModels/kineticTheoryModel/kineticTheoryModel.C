/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "kineticTheoryModel.H"
#include "mathematicalConstants.H"
#include "twoPhaseSystem.H"
#include "fvOptions.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RASModels::kineticTheoryModel::kineticTheoryModel
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& phase,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity
    <
        RASModel<EddyDiffusivity<phaseCompressibleTurbulenceModel>>
    >
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        phase,
        propertiesName
    ),

    phase_(phase),

    viscosityModel_
    (
        kineticTheoryModels::viscosityModel::New
        (
            coeffDict_
        )
    ),
    conductivityModel_
    (
        kineticTheoryModels::conductivityModel::New
        (
            coeffDict_
        )
    ),
    radialModel_
    (
        kineticTheoryModels::radialModel::New
        (
            coeffDict_
        )
    ),
    granularPressureModel_
    (
        kineticTheoryModels::granularPressureModel::New
        (
            coeffDict_
        )
    ),
    frictionalStressModel_
    (
        kineticTheoryModels::frictionalStressModel::New
        (
            coeffDict_
        )
    ),

    equilibrium_(coeffDict_.lookup("equilibrium")),
    KochSanganiDrag_(coeffDict_.lookupOrDefault<Switch>("KochSanganiDrag", true)),
    e_("e", dimless, coeffDict_),
    alphaMax_("alphaMax", dimless, coeffDict_),
    alphaMinFriction_
    (
        "alphaMinFriction",
        dimless,
        coeffDict_
    ),
    residualAlpha_
    (
        "residualAlpha",
        dimless,
        coeffDict_
    ),

    maxNut_
    (
        "maxNut",
        dimensionSet(0,2,-1,0,0),
        coeffDict_.lookupOrDefault<scalar>("maxNut",1000)
    ),

    Theta_
    (
        IOobject
        (
            IOobject::groupName("Theta", phase.name()),
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),

    lambda_
    (
        IOobject
        (
            IOobject::groupName("lambda", phase.name()),
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("zero", dimensionSet(0, 2, -1, 0, 0), 0.0)
    ),

    gs0_
    (
        IOobject
        (
            IOobject::groupName("gs0", phase.name()),
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.0)
    ),

    kappa_
    (
        IOobject
        (
            IOobject::groupName("kappa", phase.name()),
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("zero", dimensionSet(1, -1, -1, 0, 0), 0.0)
    ),

    alphaturb_
    (
        IOobject
        (
            IOobject::groupName("alphaturb", phase.name()),
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("zero", dimensionSet(1, -1, -1, 0, 0), 0.0)
    ),


    nuFric_
    (
        IOobject
        (
            IOobject::groupName("nuFric", phase.name()),
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("zero", dimensionSet(0, 2, -1, 0, 0), 0.0)
    ),

    pf_
    (
        IOobject
        (
            IOobject::groupName("pf", phase.name()),
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("zero", dimensionSet(1, -1, -2, 0, 0), 0.0)//,
        //zeroGradientFvPatchField<scalar>::typeName
     )
{
    if (type == typeName)
    {
        printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RASModels::kineticTheoryModel::~kineticTheoryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::RASModels::kineticTheoryModel::read()
{
    if
    (
        eddyViscosity
        <
            RASModel<EddyDiffusivity<phaseCompressibleTurbulenceModel>>
        >::read()
    )
    {
        coeffDict().lookup("equilibrium") >> equilibrium_;
        coeffDict().lookup("KochSanganiDrag") >> KochSanganiDrag_;
        e_.readIfPresent(coeffDict());
        alphaMax_.readIfPresent(coeffDict());
        alphaMinFriction_.readIfPresent(coeffDict());

        viscosityModel_->read();
        conductivityModel_->read();
        radialModel_->read();
        granularPressureModel_->read();
        frictionalStressModel_->read();

        return true;
    }
    else
    {
        return false;
    }
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::kineticTheoryModel::k() const
{
    return Theta_;
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::kineticTheoryModel::epsilon() const
{
    NotImplemented;
    return nut_;
}


Foam::tmp<Foam::volSymmTensorField>
Foam::RASModels::kineticTheoryModel::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                IOobject::groupName("R", U_.group()),
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
          - (nut_ + nuFric_)*dev(twoSymm(fvc::grad(U_)))
          - (lambda_*fvc::div(phi_))*symmTensor::I
        )
    );
}

void Foam::RASModels::kineticTheoryModel::boundGradU
(
    volTensorField& R
) const
{
    scalar sMin = -1.0e3;
    scalar sMax =  1.0e3;

    R.max
    (
        dimensionedTensor
        (
            "zero",
            R.dimensions(),
            tensor
            (
                  sMin, sMin, sMin,
                  sMin, sMin, sMin,
                  sMin, sMin, sMin
            )
        )
    );
    
    R.min
    (
        dimensionedTensor
        (
            "zero",
            R.dimensions(),
            tensor
            (
                  sMax, sMax, sMax,
                  sMax, sMax, sMax,
                  sMax, sMax, sMax
            )
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::kineticTheoryModel::pPrime() const
{
    const volScalarField& rho = phase_.rho();
    tmp<volScalarField> tda(phase_.d());
    const volScalarField& da = tda();
    // Get strain rate tensor for frictional pressure models
    volTensorField gradU(fvc::grad(phase_.U()));
    boundGradU(gradU);
    volSymmTensorField D(symm(gradU));
    
    tmp<volScalarField> tpPrime
    (
        Theta_
       *granularPressureModel_->granularPressureCoeffPrime
        (
            alpha_,
            radialModel_->g0(alpha_, alphaMinFriction_, alphaMax_),
            radialModel_->g0prime(alpha_, alphaMinFriction_, alphaMax_),
            rho,
            da,
            e_
        )
      + pos(alpha_ -alphaMinFriction_)
       *frictionalStressModel_->frictionalPressurePrime
        (
            phase_,
            alphaMinFriction_,
            alphaMax_,
            da,
            rho,
            D
        )
    );

    volScalarField::Boundary& bpPrime =
        tpPrime.ref().boundaryFieldRef();

    forAll(bpPrime, patchi)
    {
        if (!bpPrime[patchi].coupled())
        {
            bpPrime[patchi] = 0;
        }
    }

    return tpPrime;
}

Foam::tmp<Foam::surfaceScalarField>
Foam::RASModels::kineticTheoryModel::pPrimef() const
{
    return fvc::interpolate(pPrime());
}

Foam::tmp<Foam::volScalarField>
Foam::RASModels::kineticTheoryModel::normalStress() const
{
    const volScalarField& rho = phase_.rho();
    tmp<volScalarField> tda(phase_.d());
    const volScalarField& da = tda();
    // Get strain rate tensor for frictional pressure models
    volTensorField gradU(fvc::grad(phase_.U()));
    boundGradU(gradU);
    volSymmTensorField D(symm(gradU));
    
    tmp<volScalarField> tNormalStress
    (
        Theta_
       *granularPressureModel_->granularPressureCoeff
        (
            alpha_,
            radialModel_->g0(alpha_, alphaMinFriction_, alphaMax_),
            rho,
            da,
            e_
        )
      + pos(alpha_ - alphaMinFriction_)
       *frictionalStressModel_->frictionalPressure
        (
            phase_,
            alphaMinFriction_,
            alphaMax_,
            da,
            rho,
            D
        )
      );

    return tNormalStress;
}


Foam::tmp<Foam::volSymmTensorField>
Foam::RASModels::kineticTheoryModel::devRhoReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                IOobject::groupName("devRhoReff", U_.group()),
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
          - (rho_*(nut_ + nuFric_))
           *dev(twoSymm(fvc::grad(U_)))
          - ((rho_*lambda_)*fvc::div(phi_))*symmTensor::I
        )
    );
}


Foam::tmp<Foam::fvVectorMatrix>
Foam::RASModels::kineticTheoryModel::divDevRhoReff
(
    volVectorField& U
) const
{
    return
    (
      - fvm::laplacian(rho_*(nut_ + nuFric_), U)
      - fvc::div
        (
            (rho_*(nut_ + nuFric_))*dev2(T(fvc::grad(U)))
          + ((rho_*lambda_)*fvc::div(phi_))
           *dimensioned<symmTensor>("I", dimless, symmTensor::I)
        )
    );
}


void Foam::RASModels::kineticTheoryModel::correct()
{
    // Local references
    const twoPhaseSystem& fluid = refCast<const twoPhaseSystem>(phase_.fluid());
    volScalarField alpha(max(alpha_, residualAlpha_));
    const volScalarField& rho = phase_.rho();
    const surfaceScalarField& alphaRhoPhi = alphaRhoPhi_;
    const volVectorField& U = U_;
    const volVectorField& Uc_ = fluid.otherPhase(phase_).U();

    const scalar sqrtPi = sqrt(constant::mathematical::pi);
    dimensionedScalar ThetaSmall("ThetaSmall", Theta_.dimensions(), 1.0e-8);
    dimensionedScalar ThetaSmallSqrt(sqrt(ThetaSmall));

    tmp<volScalarField> tda(phase_.d());
    const volScalarField& da = tda();

    volTensorField gradU(fvc::grad(U_));
    boundGradU(gradU);
    
    volSymmTensorField D(symm(gradU));
    
    volScalarField cellVolume
    (
        IOobject
        (
            "cellVolume",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("one", dimLength*dimLength*dimLength, 1)
    );
    cellVolume.ref() = mesh_.V();

    // Calculating the radial distribution function
    gs0_ = radialModel_->g0(alpha, alphaMinFriction_, alphaMax_);
    
    volScalarField trD(fvc::div(phi_));
    
    // Drag
    const volScalarField& beta = mesh_.lookupObject<volScalarField>("Kd");
    
    volScalarField ThetaSqrt("sqrtTheta", sqrt(Theta_));
    
    volScalarField PsCoeff
    (
        granularPressureModel_->granularPressureCoeff
        (
            alpha,
            gs0_,
            rho,
            da,
            e_
        )
    );
    
    Info << "Kinetic Theory:" << endl;
    
    if (!equilibrium_) {
        // Particle viscosity (Table 3.2, p.47)
        nut_ = viscosityModel_->nu(alpha, Theta_, gs0_, rho, da, e_);

        // Bulk viscosity  p. 45 (Lun et al. 1984).
        lambda_ = (4.0/3.0)*sqr(alpha)*da*gs0_*(1.0 + e_)*ThetaSqrt/sqrtPi;

        // Dissipation (Eq. 3.24, p.50)
        volScalarField gammaCoeff
        (
            "gammaCoeff",
            12.0
           *(1.0 - sqr(e_))
           *max(sqr(alpha), sqr(residualAlpha_))
           *rho
           *gs0_
           *ThetaSqrt
           /(da*sqrtPi)
        );

        // Gidaspow (1994)
        // volScalarField J1("J1", 3.0*beta);
        volScalarField J1
        (
            "J1",
            3.0*beta
        );
        volScalarField J2
        (
            "J2",
           0.25*sqr(beta)*da*magSqr(U - Uc_)
           /(
               rho
              *max(alpha,residualAlpha_)
              *radialModel_->g0(alpha, 0.999*alphaMax_, alphaMax_)
              *sqrtPi
              *(ThetaSqrt*Theta_ + ThetaSmallSqrt*ThetaSmall)
            )
        );
        if (KochSanganiDrag_) {
            // Koch & Sangani (1999): drag dissipation
            volScalarField Rdiss
            (
                1.0
              + 3.0*sqrt(0.5*alpha)
              + 2.11*alpha*log(alpha)
              + 11.26*alpha*(1.0 - 5.1*alpha + 16.57*sqr(alpha) - 21.77*pow3(alpha))
              - alpha*gs0_*log(0.01)
            );
            J1 *= Rdiss;
            // Koch & Sangani (1999): drag production
            volScalarField Rd
            (
                neg(alpha - 0.4)
               *(1.0 + 3.0*sqrt(0.5*alpha) + 2.11*alpha*log(alpha)+17.14*alpha)
               /(1.0 + 0.681*alpha - 8.48*sqr(alpha) + 8.16*pow3(alpha))
              + pos(alpha - 0.4)
               *(10.0*alpha/(pow3((1.0 - alpha))) + 0.7)
            );
            volScalarField Psi
            (
                sqr(Rd)
               /(1.0 + 2.5*sqrt(alpha) + 5.9*alpha)
            );
            J2 *= Psi;            
        }

        // 'thermal' conductivity (Table 3.3, p. 49)
        kappa_ = conductivityModel_->kappa(alpha, Theta_, gs0_, rho, da, e_);

        fv::options& fvOptions(fv::options::New(mesh_));

        // Construct the granular temperature equation (Eq. 3.20, p. 44)
        fvScalarMatrix ThetaEqn
        (
            1.5
           *(
                fvm::ddt(alpha, rho, Theta_)
              + fvm::div(alphaRhoPhi, Theta_)
              + fvm::SuSp(-(fvc::ddt(alpha, rho) + fvc::div(alphaRhoPhi)), Theta_)
            )
          - fvm::laplacian(kappa_, Theta_, "laplacian(kappa,Theta)")
         ==
          - fvm::SuSp(PsCoeff*trD, Theta_)
          + rho
           *(
                lambda_*sqr(trD)
              + 2.0*nut_*(dev(D)&&gradU)
            )
          + fvm::Sp(-gammaCoeff, Theta_)
          + fvm::Sp(-J1, Theta_)
          + fvm::Sp(J2, Theta_)
          + fvOptions(alpha, rho, Theta_)
        );

        ThetaEqn.relax();
        fvOptions.constrain(ThetaEqn);
        ThetaEqn.solve();
        fvOptions.correct(Theta_);
    } else {
        // Equilibrium => dissipation == production
        // Eq. 4.14, p.82
        volScalarField K1("K1", 2.0*(1.0 + e_)*rho*gs0_);
        volScalarField K3
        (
            "K3",
            0.5*da*rho*
            (
                (sqrtPi/(3.0*(3.0 - e_)))
               *(1.0 + 0.4*(1.0 + e_)*(3.0*e_ - 1.0)*alpha*gs0_)
               +1.6*alpha*gs0_*(1.0 + e_)/sqrtPi
            )
        );

        volScalarField K2
        (
            "K2",
            4.0*da*rho*(1.0 + e_)*alpha*gs0_/(3.0*sqrtPi) - 2.0*K3/3.0
        );

        volScalarField K4("K4", 12.0*(1.0 - sqr(e_))*rho*gs0_/(da*sqrtPi));

        volScalarField trD
        (
            "trD",
            alpha/(alpha + residualAlpha_)
           *fvc::div(phi_)
        );
        volScalarField tr2D("tr2D", sqr(trD));
        volScalarField trD2("trD2", tr(D & D));

        volScalarField t1("t1", K1*alpha + rho);
        volScalarField l1("l1", -t1*trD);
        volScalarField l2("l2", sqr(t1)*tr2D);
        volScalarField l3
        (
            "l3",
            4.0
           *K4
           *alpha
           *(2.0*K3*trD2 + K2*tr2D)
        );

        Theta_ = sqr
        (
            (l1 + sqrt(l2 + l3))
           /(2.0*max(alpha, residualAlpha_)*K4)
        );

        kappa_ = conductivityModel_->kappa(alpha, Theta_, gs0_, rho, da, e_);
    }

    Theta_.max(ThetaSmall.value());
    Theta_.min(100);
    
    Info << "    max(Theta) = " << max(Theta_).value() << endl;
    {
        // particle viscosity (Table 3.2, p.47)
        nut_ = viscosityModel_->nu(alpha, Theta_, gs0_, rho, da, e_);

        ThetaSqrt = sqrt(Theta_);

        // Bulk viscosity  p. 45 (Lun et al. 1984).
        lambda_ = (4.0/3.0)*sqr(alpha)*da*gs0_*(1.0 + e_)*ThetaSqrt/sqrtPi;

        // Frictional pressure
        pf_ = frictionalStressModel_->frictionalPressure
        (
            phase_,
            alphaMinFriction_,
            alphaMax_,
            da,
            rho,
            D
        );
        //pf_.correctBoundaryConditions();

        nuFric_ = frictionalStressModel_->nu
        (
            phase_,
            alphaMinFriction_,
            alphaMax_,
            pf_/rho,
            da,
            D
        );
        
        // Limit viscosity and add frictional viscosity
        nut_.min(maxNut_);
        nuFric_.min(maxNut_);

        alphaturb_ = rho*da*pow(constant::mathematical::pi,1.5)*sqrt(Theta_)/(32*(16 - 7*alpha)/(16*(1 - alpha)*(1 - alpha)));
        
        
        Info << "    max(nut) = " << max(nut_).value() << nl
             << "    max(nutFric) = " << max(nuFric_).value() << endl;
        

    }

    if (debug)
    {
        Info<< typeName << ':' << nl
            << "    max(Theta) = " << max(Theta_).value() << nl
            << "    max(nut) = " << max(nut_).value() << endl;
    }
}


// ************************************************************************* //

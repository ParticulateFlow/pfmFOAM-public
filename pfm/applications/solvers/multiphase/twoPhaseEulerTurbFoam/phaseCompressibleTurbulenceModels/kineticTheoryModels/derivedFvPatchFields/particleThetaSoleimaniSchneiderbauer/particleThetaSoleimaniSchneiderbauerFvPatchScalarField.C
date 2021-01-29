/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 
     \\/     M anipulation  | VTT Technical Research Centre of Finland Ltd
-------------------------------------------------------------------------------
License
    This file is a derived work of OpenFOAM.

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

#include "particleThetaSoleimaniSchneiderbauerFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

#include "twoPhaseSystem.H"
#include "compressibleTurbulenceModel.H"
#include "ThermalDiffusivity.H"
#include "PhaseCompressibleTurbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

particleThetaSoleimaniSchneiderbauerFvPatchScalarField::
particleThetaSoleimaniSchneiderbauerFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<double, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    restitutionCoefficient_(0),
    tangentialRestitutionCoeff_(0),
    muF_(0),
    sigma_(0),
    residualAlpha_(1e-8),
    residualKappa_(1e-15),
    residualTheta_(1e-8),
    residualU_(1e-8),
    mu0_(0),
    eta_(0)  
{}


particleThetaSoleimaniSchneiderbauerFvPatchScalarField::
particleThetaSoleimaniSchneiderbauerFvPatchScalarField
(
    const particleThetaSoleimaniSchneiderbauerFvPatchScalarField& ptsssf,
    const fvPatch& p,
    const DimensionedField<double, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptsssf, p, iF, mapper),
    restitutionCoefficient_(ptsssf.restitutionCoefficient_),
    tangentialRestitutionCoeff_(ptsssf.tangentialRestitutionCoeff_),
    muF_(ptsssf.muF_),
    sigma_(ptsssf.sigma_),
    residualAlpha_(ptsssf.residualAlpha_),
    residualKappa_(ptsssf.residualKappa_),
    residualTheta_(ptsssf.residualTheta_),
    residualU_(ptsssf.residualU_),
    mu0_(ptsssf.mu0_),
    eta_(ptsssf.eta_)    
{}


particleThetaSoleimaniSchneiderbauerFvPatchScalarField::
particleThetaSoleimaniSchneiderbauerFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<double, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    restitutionCoefficient_(readScalar(dict.lookup("restitutionCoefficient"))),
    tangentialRestitutionCoeff_(readScalar(dict.lookup("tangentialRestitutionCoeff"))),
    muF_(readScalar(dict.lookup("muF"))),
    sigma_(readScalar(dict.lookup("sigma"))),
    residualAlpha_(dict.lookupOrDefault<scalar>("residualAlpha",1e-8)),
    residualKappa_(dict.lookupOrDefault<scalar>("residualKappa",1e-12)),
    residualTheta_(dict.lookupOrDefault<scalar>("residualTheta",1e-8)),
    residualU_(dict.lookupOrDefault<scalar>("residualU",1e-8)),
    mu0_(7.0/2.0*(1.0 + restitutionCoefficient_)/(1.0 + tangentialRestitutionCoeff_)*muF_),
    eta_(1.0/2.0*(1.0 + restitutionCoefficient_))   
{
    if (dict.found("value") && dict.found("gradient"))
    {
        fvPatchField<scalar>::operator=(Field<scalar>("value", dict, p.size()));
        gradient() = Field<scalar>("gradient", dict, p.size());
    }
    else
    {
        // Still reading so cannot yet evaluate. Make up a value.
        fvPatchField<scalar>::operator=(patchInternalField());
        gradient() = 0.0;
    }
}


particleThetaSoleimaniSchneiderbauerFvPatchScalarField::
particleThetaSoleimaniSchneiderbauerFvPatchScalarField
(
    const particleThetaSoleimaniSchneiderbauerFvPatchScalarField& ptsssf,
    const DimensionedField<double, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(ptsssf, iF),
    restitutionCoefficient_(ptsssf.restitutionCoefficient_),
    tangentialRestitutionCoeff_(ptsssf.tangentialRestitutionCoeff_),
    muF_(ptsssf.muF_),
    sigma_(ptsssf.sigma_),
    residualAlpha_(ptsssf.residualAlpha_),
    residualKappa_(ptsssf.residualKappa_),
    residualTheta_(ptsssf.residualTheta_),
    residualU_(ptsssf.residualU_),
    mu0_(ptsssf.mu0_),
    eta_(ptsssf.eta_)   
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void particleThetaSoleimaniSchneiderbauerFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    scalarField::autoMap(m);
}


void particleThetaSoleimaniSchneiderbauerFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchScalarField::rmap(ptf, addr);
}


// Update the coefficients associated with the patch field
void particleThetaSoleimaniSchneiderbauerFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    if ((restitutionCoefficient_ < 0) || (restitutionCoefficient_ > 1))
    {
        FatalErrorInFunction
            << "The specularity coefficient has to be between 0 and 1"
            << abort(FatalError);
    }

    const label patchi = patch().index();

    // lookup the fluid model and the phase
    const twoPhaseSystem& fluid = db().lookupObject<twoPhaseSystem>
    (
        "phaseProperties"
    );

    const phaseModel& granular
    (
        fluid.phase1().name() == internalField().group()
      ? fluid.phase1()
      : fluid.phase2()
    );
    // lookup all the fields on this patch
    const fvPatchScalarField& alpha
    (
        patch().lookupPatchField<volScalarField, scalar>
        (
            granular.volScalarField::name()
        )
    );
    const scalarField alphap(max(alpha, residualAlpha_));

    // Calculate velocity components:
    const tmp<volVectorField> tU = granular.U();
    const volVectorField& U = tU();
    const fvPatchVectorField& Up = U.boundaryField()[patchi];
    const vectorField Uc(Up.patchInternalField());
    const scalarField magUc(mag(Uc));

    // Wall normal cell velocity
    const scalarField Un(-(this->patch().nf() & Uc));

    // Tangential cell velocity
    const scalarField Utc(mag(Uc + this->patch().nf()*Un));

    const scalarField rhop
    (
        patch().lookupPatchField<volScalarField, scalar>
        (
            IOobject::groupName("thermo:rho", granular.name())
        )
    );
    const scalarField Thetap(max(*this, residualTheta_));

    // Lookup kinetic theory properties not available trhough the 
    // turbulenceModel interface

    const fvPatchScalarField& gs0p
    (
        patch().lookupPatchField<volScalarField, scalar>
        (
            IOobject::groupName("gs0", granular.name())
        )
    );

    const fvPatchScalarField& kappap
    (
        patch().lookupPatchField<volScalarField, scalar>
        (
            IOobject::groupName("kappa", granular.name())
        )
    );

    // lookup the packed volume fraction
    dimensionedScalar alphaMax
    (
        "alphaMax",
        dimless,
        db()
       .lookupObject<IOdictionary>
        (
            IOobject::groupName("turbulenceProperties", granular.name())
        )
       .subDict("RAS")
       .subDict("kineticTheoryCoeffs")
       .lookup("alphaMax")
    );

    // Impact angle
    const scalarField Imp
    (
        radToDeg
        (
            atan
            (
                max(sqrt(1.5*Thetap)/2 + Un, scalar(0.0))
               /max(Utc, residualU_)
            )
        )
    );

    // Virtual wall angle
    const scalarField gamma
    (
        degToRad
        (
            2*0.398942*sigma_
           *(
                exp(-0.5*sqr(Imp)/(sqr(sigma_)))
              - exp(-4050/(sqr(sigma_)))
            )
           /(erf(63.6396/sigma_) + erf(0.707107*Imp/sigma_))
        )
    );

    const scalarField u
    (
        Utc*cos(gamma)
      - Un*sin(gamma)
    );

    const scalarField v
    (
       - Utc*sin(gamma)
       + Un*cos(gamma)
    );

    const scalarField us((u + mu0_*v)/(sqrt(2.0*Thetap)*mu0_));

    const scalarField erfUs(erf(us));

    // Wall shear stress
    const scalarField tauw
    ( 
      - rhop*alphap*gs0p*eta_*muF_
       *(
            u*v*(1-erfUs)/mu0_
          + v*sqrt(2*Thetap/constant::mathematical::pi)*
            (
                exp(-sqr(v)/(2*Thetap))
              - exp(-sqr(us))
            )
          + (sqr(v)+Thetap)
           *(
                erf(v/sqrt(2*Thetap))
              - erfUs
            )
        )
    );

    // Wall normal stress
    const scalarField N
    (
        rhop*alphap*gs0p*eta_
       *(
            (sqr(v) + Thetap)*(1 - erf(v/sqrt(2*Thetap)))
          - sqrt(2*Thetap/constant::mathematical::pi)
           *v*exp(-sqr(v)/(2*Thetap))
        )
    );

    const scalarField A
    (
        2*muF_*sqr(u)*(2*eta_ - mu0_)
      + Thetap
       *(
            14.0*muF_*eta_
          - 4*mu0_*(1.0 + muF_)
          - 6.0*muF_*sqr(mu0_)*eta_
        )
    );

    const scalarField Ul(u/(sqrt(2*Thetap)*mu0_));

    const scalarField B
    (
        sqrt(Thetap)
       *(
            4.0*(eta_ - 1.0)
          + 6.0*sqr(muF_)*eta_
        )
      - sqrt(2.0*constant::mathematical::pi)*muF_*u*erf(Ul)
    );
/*
    const scalarField D01
    (
      - alphap*rhop*gs0p*eta_*sqrt(Thetap)
       /(sqr(mu0_)*sqrt(2.0*constant::mathematical::pi))
       *(
            exp(-sqr(Ul))*(muF_*A)
        )
    );
*/
    const scalarField D0
    (
      - alphap*rhop*gs0p*eta_*sqrt(Thetap)
       /(sqr(mu0_)*sqrt(2.0*constant::mathematical::pi))
       *(
            exp(-sqr(Ul))*(muF_*A)
          + sqr(mu0_)*sqrt(Thetap)*B
        )
    );

    const scalarField N0(eta_*rhop*alphap*gs0p*Thetap);
    const scalarField D(D0/N0*N);
    
    this->gradient() = -((tauw*(u) - N*(v)) + D)/(max(kappap, residualKappa_));

    // Info<< "  theta-grad: " << min(this->gradient()) << " - " << max(this->gradient()) << endl;

    fixedGradientFvPatchScalarField::updateCoeffs();
}


// Write
void particleThetaSoleimaniSchneiderbauerFvPatchScalarField::
write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeKeyword("restitutionCoefficient")
        << restitutionCoefficient_ << token::END_STATEMENT << nl;
    os.writeKeyword("tangentialRestitutionCoeff")
        << tangentialRestitutionCoeff_ << token::END_STATEMENT << nl;
    os.writeKeyword("muF")
        << muF_ << token::END_STATEMENT << nl;
    os.writeKeyword("sigma")
        << sigma_ << token::END_STATEMENT << nl;
    writeEntry("value",os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    particleThetaSoleimaniSchneiderbauerFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

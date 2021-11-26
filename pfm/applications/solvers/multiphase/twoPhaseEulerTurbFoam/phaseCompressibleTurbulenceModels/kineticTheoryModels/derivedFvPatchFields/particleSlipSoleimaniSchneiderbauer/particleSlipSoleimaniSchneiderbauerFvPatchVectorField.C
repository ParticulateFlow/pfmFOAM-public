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

#include "particleSlipSoleimaniSchneiderbauerFvPatchVectorField.H"
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

particleSlipSoleimaniSchneiderbauerFvPatchVectorField::
particleSlipSoleimaniSchneiderbauerFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<Vector<double>, volMesh>& iF
)
:
    partialSlipFvPatchVectorField(p, iF),
    restitutionCoefficient_(0),
    tangentialRestitutionCoeff_(0),
    muF_(0),
    sigma_(0),
    residualAlpha_(1.0e-8),
    residualTheta_(1.0e-8),
    residualU_(1.0e-8),
    mu0_(0)
 {}


particleSlipSoleimaniSchneiderbauerFvPatchVectorField::
particleSlipSoleimaniSchneiderbauerFvPatchVectorField
(
    const particleSlipSoleimaniSchneiderbauerFvPatchVectorField& psssvf,
    const fvPatch& p,
    const DimensionedField<Vector<double>, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    partialSlipFvPatchVectorField(psssvf, p, iF, mapper),
    restitutionCoefficient_(psssvf.restitutionCoefficient_),
    tangentialRestitutionCoeff_(psssvf.tangentialRestitutionCoeff_),
    muF_(psssvf.muF_),
    sigma_(psssvf.sigma_),
    residualAlpha_(psssvf.residualAlpha_),
    residualTheta_(psssvf.residualTheta_),
    residualU_(psssvf.residualU_),
    mu0_(psssvf.mu0_)
   
{}


particleSlipSoleimaniSchneiderbauerFvPatchVectorField::
particleSlipSoleimaniSchneiderbauerFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<Vector<double>, volMesh>& iF,
    const dictionary& dict
)
:
    partialSlipFvPatchVectorField(p, iF),
    restitutionCoefficient_(readScalar(dict.lookup("restitutionCoefficient"))),
    tangentialRestitutionCoeff_(readScalar(dict.lookup("tangentialRestitutionCoeff"))),
    muF_(readScalar(dict.lookup("muF"))),
    sigma_(readScalar(dict.lookup("sigma"))),
    residualAlpha_(dict.lookupOrDefault<scalar>("residualAlpha",1.0e-8)),
    residualTheta_(dict.lookupOrDefault<scalar>("residualTheta",1.0e-8)),
    residualU_(dict.lookupOrDefault<scalar>("residualU",1.0e-8)),
    mu0_(7.0/2.0*(1.0 + restitutionCoefficient_)/(1.0 + tangentialRestitutionCoeff_)*muF_)
   
{
    if (dict.found("value"))
    {
        fvPatchField<vector>::operator=
        (
            vectorField("value", dict, p.size())
        );
    }
    else
    {
        partialSlipFvPatchVectorField::evaluate();
    }
}


particleSlipSoleimaniSchneiderbauerFvPatchVectorField::
particleSlipSoleimaniSchneiderbauerFvPatchVectorField
(
    const particleSlipSoleimaniSchneiderbauerFvPatchVectorField& psssvf,
    const DimensionedField<Vector<double>, volMesh>& iF
)
:
    partialSlipFvPatchVectorField(psssvf, iF),
    restitutionCoefficient_(psssvf.restitutionCoefficient_),
    tangentialRestitutionCoeff_(psssvf.tangentialRestitutionCoeff_),
    muF_(psssvf.muF_),
    sigma_(psssvf.sigma_),
    residualAlpha_(psssvf.residualAlpha_),
    residualTheta_(psssvf.residualTheta_),
    residualU_(psssvf.residualU_),
    mu0_(psssvf.mu0_)
  
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
// Update the coefficients associated with the patch field
void particleSlipSoleimaniSchneiderbauerFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
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

    // Wall normal cell velocity
    const scalarField Un(-(this->patch().nf() & Uc));

    // Tangential cell velocity
    const scalarField Utc(mag(Uc + this->patch().nf()*Un));
    const scalarField magUc(mag(Utc)+scalar(1.0e-8));
    
    const scalarField rhop
    (
        patch().lookupPatchField<volScalarField, scalar>
        (
            IOobject::groupName("thermo:rho", granular.name())
        )
    );
    
    const scalarField pfW
    (
        patch().lookupPatchField<volScalarField, scalar>
        (
            IOobject::groupName("pf", granular.name())
         )
    );
    const scalarField nu
    (
        patch().lookupPatchField<volScalarField, scalar>
        (
            IOobject::groupName("nut", granular.name())
        )
    );
    
    const scalarField nuF
    (
        patch().lookupPatchField<volScalarField, scalar>
        (
            IOobject::groupName("nuFric", granular.name())
        )
    );
    const scalarField muP = nu*rhop;
    const scalarField muFP = nuF*rhop;

    word ThetaName(IOobject::groupName("Theta", granular.name()));

    const fvPatchScalarField& Theta
    (
        db().foundObject<volScalarField>(ThetaName)
      ? patch().lookupPatchField<volScalarField, scalar>(ThetaName)
      : alpha
    );

    const scalarField Thetap(max(Theta, residualTheta_));

    const fvPatchScalarField& gs0p
    (
        patch().lookupPatchField<volScalarField, scalar>
        (
            IOobject::groupName("gs0", granular.name())
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
                max(sqrt(1.5*Thetap)/2.0 + Un, scalar(0.0))
               /max(Utc, residualU_)
            )
        )
    );

    // Virtual wall angle
    const scalarField gamma
    (
        degToRad
        (
            2.0*0.398942*sigma_
           *(
                exp(-0.5*sqr(Imp)/(sqr(sigma_ + scalar(1.0e-8))))
              - exp(-4050/(sqr(sigma_ + scalar(1.0e-8))))
            )
           /(erf(63.6396/(sigma_ + scalar(1.0e-8))) + erf(0.707107*Imp/(sigma_ + scalar(1.0e-8))))
        )
    );

    // u and v should have the magnitude of the particle
    // velocity at the wall
    const scalarField u
    (
        (
            Utc*cos(gamma)
          - Un*sin(gamma)
        )
//       *mag(Up)
//       /magUc
    );

    const scalarField v
    (
        (
          - Utc*sin(gamma)
          + Un*cos(gamma)
        )
//       *mag(Up)
//       /magUc
    ); 

    const scalarField us
    (
        (u + mu0_*v)
       /(sqrt(2.0*Thetap)*mu0_)
    );

    const scalarField erfUs(erf(us));

    // Wall shear stress
    const scalarField tauw
    (
        0.5*rhop*alphap*gs0p*(1.0 + restitutionCoefficient_)*muF_
       *(
             u*v*(1.0-erfUs)/mu0_
          + v*sqrt(2.0*Thetap/constant::mathematical::pi)*
            (
                exp(-sqr(v)/(2.0*Thetap))
              - exp(-sqr(us))
            )
          + (sqr(v)+Thetap)
           *(
                erf(v/sqrt(2.0*Thetap))
              - erfUs
            )
        ) 
    );

    /*=======================================================================*\

       The partial slip BC in OpenFOAM is implemented as

       (1) Up = (1 - valueFraction)*Uc

       Here Up is velocity on the patch and Uc is velocity in the 
       first cell centre.

       Wall shear stress is: 

       (2) tauw = muEffp*(Up-Uc)*deltaCoeff

       To set the tauw calculated above, the above equation is solved for Up:

       (3) Up = tauw/(muEffp*deltaCoeff)+Uc

       Inserting this into (1) and solving for valueFraction:

       (4) valueFraction = -tauw/(muEffp*Uc*deltaCoeff)
     
       However, the value fraction should be 0-1 and thus the maximum shear
       that can be bescribed with the present partialSlip based formulation is

       (5) tauw = -muEffp*Uc*deltaCoeff

       If the above tauw correlation predicts a higher wall shear stress 
       than (5) it is neglected and the velocity at the wall is set to zero.   

    \*=======================================================================*/    

    this->valueFraction() = 
        max
        (
            min
            (
                (-tauw + 0.0*pfW*muF_)/max(Utc*(muP+0.0*muFP)*patch().deltaCoeffs(),scalar(1.0e-8)),
                scalar(1.0)
            ),
            scalar(0)
        );

    Info<< "  tauw: " << min(tauw) << " - " << max(tauw) << endl;
    Info<< "  pfW: "  << min(pfW) << " - " << max(pfW) << endl;
    Info<< "  valueFraction(): " << min(this->valueFraction()) 
        << " - " << max(this->valueFraction()) << endl;
        
    partialSlipFvPatchVectorField::updateCoeffs();       
}

// Write
void particleSlipSoleimaniSchneiderbauerFvPatchVectorField::
write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    os.writeKeyword("restitutionCoefficient")
        << restitutionCoefficient_ << token::END_STATEMENT << nl;
    os.writeKeyword("tangentialRestitutionCoeff")
        << tangentialRestitutionCoeff_ << token::END_STATEMENT << nl;
    os.writeKeyword("muF")
        << muF_ << token::END_STATEMENT << nl;
    os.writeKeyword("sigma")
        << sigma_ << token::END_STATEMENT << nl;
    writeEntry("value", os);

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    particleSlipSoleimaniSchneiderbauerFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

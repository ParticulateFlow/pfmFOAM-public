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

#include "partialSlipFrictionFvPatchVectorField.H"
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

partialSlipFrictionFvPatchVectorField::
partialSlipFrictionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<Vector<double>, volMesh>& iF
)
:
    partialSlipFvPatchVectorField(p, iF),
    muW_(0),
    residualAlpha_(1e-8),
    residualU_(1e-8)
 {}


partialSlipFrictionFvPatchVectorField::
partialSlipFrictionFvPatchVectorField
(
    const partialSlipFrictionFvPatchVectorField& psssvf,
    const fvPatch& p,
    const DimensionedField<Vector<double>, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    partialSlipFvPatchVectorField(psssvf, p, iF, mapper),
    muW_(psssvf.muW_),
    residualAlpha_(psssvf.residualAlpha_),
    residualU_(psssvf.residualU_)
{}


partialSlipFrictionFvPatchVectorField::
partialSlipFrictionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<Vector<double>, volMesh>& iF,
    const dictionary& dict
)
:
    partialSlipFvPatchVectorField(p, iF),
    muW_(readScalar(dict.lookup("muW"))),
    residualAlpha_(dict.lookupOrDefault<scalar>("residualAlpha",1e-8)),
    residualU_(dict.lookupOrDefault<scalar>("residualU",1e-8))
   
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


partialSlipFrictionFvPatchVectorField::
partialSlipFrictionFvPatchVectorField
(
    const partialSlipFrictionFvPatchVectorField& psssvf,
    const DimensionedField<Vector<double>, volMesh>& iF
)
:
    partialSlipFvPatchVectorField(psssvf, iF),
    muW_(psssvf.muW_),
    residualAlpha_(psssvf.residualAlpha_),
    residualU_(psssvf.residualU_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Update the coefficients associated with the patch field
void partialSlipFrictionFvPatchVectorField::updateCoeffs()
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
    const scalarField Un(this->patch().nf() & Uc);

    // Tangential cell velocity
    const scalarField Utc(mag(Uc - this->patch().nf()*Un));

    const vectorField kp
    (
        patch().lookupPatchField<volVectorField, vector>
        (
            IOobject::groupName("k", granular.name())
        )
    );
    const scalarField kpn(mag(kp - (this->patch().nf() & kp)*(this->patch().nf())));
    
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
    const scalarField muFP = nuF*rhop;
    
    const scalarField c
    (
        alphap
       *sqrt(2.0*kpn)
       *muW_
       /max(6.0*nu, VSMALL)
    );
    const scalarField tauw = -c*nu*rhop;

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

    this->valueFraction() = c/(c + patch().deltaCoeffs());
    /*    max
        (
            min
            (
                c/(c + patch().deltaCoeffs())
              + pfW*muW_
               /max(Utc*(muFP)*patch().deltaCoeffs(),scalar(1e-10)),
                scalar(1)
            ),
            scalar(0)
        );
     */
    Info<< "  tauW: "  << min(tauw) << " - " << max(tauw) << endl;
    Info<< "  pfW: "   << min(pfW) << " - " << max(pfW) << endl;
    Info<< "  valueFraction(): " << min(this->valueFraction()) 
        << " - " << max(this->valueFraction()) << endl;
        
    partialSlipFvPatchVectorField::updateCoeffs();       
}

// Write
void partialSlipFrictionFvPatchVectorField::
write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    os.writeKeyword("muW")
        << muW_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    partialSlipFrictionFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

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

#include "particleKFvPatchVectorField.H"
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

particleKFvPatchVectorField::
particleKFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<Vector<double>, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(p, iF),
    muW_(0),
    CepsW_(0),
    sigma_(1),
    residualAlpha_(1e-7)
{}


particleKFvPatchVectorField::
particleKFvPatchVectorField
(
    const particleKFvPatchVectorField& ptsssf,
    const fvPatch& p,
    const DimensionedField<Vector<double>, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchVectorField(ptsssf, p, iF, mapper),
    muW_(ptsssf.muW_),
    CepsW_(ptsssf.CepsW_),
    sigma_(ptsssf.sigma_),
    residualAlpha_(ptsssf.residualAlpha_)
{}


particleKFvPatchVectorField::
particleKFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<Vector<double>, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchVectorField(p, iF),
    muW_(readScalar(dict.lookup("muW"))),
    CepsW_(readScalar(dict.lookup("CepsW"))),
    sigma_(dict.lookupOrDefault<scalar>("sigma",1)),
    residualAlpha_(dict.lookupOrDefault<scalar>("residualAlpha",1e-7))
{
    if (dict.found("value") && dict.found("gradient"))
    {
        fvPatchField<vector>::operator=(Field<vector>("value", dict, p.size()));
        gradient() = Field<vector>("gradient", dict, p.size());
    }
    else
    {
        // Still reading so cannot yet evaluate. Make up a value.
        fvPatchField<vector>::operator=(patchInternalField());
        gradient() = vector(0.0,0.0,0.0);
    }
}


particleKFvPatchVectorField::
particleKFvPatchVectorField
(
    const particleKFvPatchVectorField& ptsssf,
    const DimensionedField<Vector<double>, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(ptsssf, iF),
    muW_(ptsssf.muW_),
    CepsW_(ptsssf.CepsW_),
    sigma_(ptsssf.sigma_),
    residualAlpha_(ptsssf.residualAlpha_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void particleKFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    vectorField::autoMap(m);
}


void particleKFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchVectorField::rmap(ptf, addr);
}


// Update the coefficients associated with the patch field
void particleKFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    Vector<double> eX(1,0,0);
    Vector<double> eY(0,1,0);
    Vector<double> eZ(0,0,1);

    
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
    const vectorField Un(-(this->patch().nf() & Uc)*this->patch().nf());

    // Tangential cell velocity
    const vectorField Utp(Uc + Un);

    const vectorField kp
    (
        patch().lookupPatchField<volVectorField, vector>
        (
            IOobject::groupName("k", granular.name())
        )
    );
    const scalarField kpn(mag(this->patch().nf() & kp));
    
    const scalarField rhop
    (
        patch().lookupPatchField<volScalarField, scalar>
        (
            IOobject::groupName("thermo:rho", granular.name())
        )
    );

    const scalarField nu
    (
        patch().lookupPatchField<volScalarField, scalar>
        (
            IOobject::groupName("nut", granular.name())
        )
    );
    
    const scalarField kappap
    (
        rhop*nu/sigma_
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
    
    const vectorField tauwUD
    (
        constant::mathematical::pi
       *alphap
       *rhop
       *(
            muW_
           *sqrt(2.0*kpn)
           *(
                (sqr(Utp&eX))*eX
              + (sqr(Utp&eY))*eY
              + (sqr(Utp&eZ))*eZ
            )
           /6.0
          - CepsW_
           *sqrt(2.0*mag(kp))
           *(
                mag(kp&eX)*eX
              + mag(kp&eY)*eY
              + mag(kp&eZ)*eZ
            )
           /4.0
         )
        /alphaMax.value()
    );
    
    this->gradient() = tauwUD/(max(kappap, SMALL));
    /*
    Info << "particleKBC: "
         << "    k-grad: "
         << min(this->gradient())
         << " - "
         << max(this->gradient())
         << endl;
    Info << "particleKBC: "
         << "    k: "
         << min(kpn)
         << " - "
         << max(kpn)
         << endl;
    */
    fixedGradientFvPatchVectorField::updateCoeffs();
}


// Write
void particleKFvPatchVectorField::
write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    os.writeKeyword("muW")
        << muW_ << token::END_STATEMENT << nl;
    os.writeKeyword("CepsW")
        << CepsW_ << token::END_STATEMENT << nl;
    os.writeKeyword("sigma")
        << sigma_ << token::END_STATEMENT << nl;
    writeEntry("value",os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    particleKFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

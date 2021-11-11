/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenFOAM Foundation
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

 (c) Simon Schneiderbauer 2021
     Christian Doppler Laboratory for Multiscale Modeling of Multiphase Processes
     Johannes Kepler University, Linz, Austria

 Description
     Combination of Schneiderbauer et al. and Chialvo et al. frictional model

     Reference:
     \verbatim
        Chialvo, S., J. Sun, S. Sundaresan. Phys. Rev. E, 2012, 85, 021305 (2012).
        Schneiderbauer, S., A. Aigner, S. Pirker. Chem. Eng. Sci., 2012, 80, 279–292.
     \endverbatim
 
 

\*---------------------------------------------------------------------------*/

#include "SchneiderbauerEtAlOldFrictionalStress.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace frictionalStressModels
{
    defineTypeNameAndDebug(SchneiderbauerEtAlOld, 0);

    addToRunTimeSelectionTable
    (
        frictionalStressModel,
        SchneiderbauerEtAlOld,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::frictionalStressModels::
SchneiderbauerEtAlOld::SchneiderbauerEtAlOld
(
    const dictionary& dict
)
:
    frictionalStressModel(dict),
    coeffDict_(dict.optionalSubDict(typeName + "Coeffs")),
    b_("b", dimless, coeffDict_),
    muSt_("muSt", dimless, coeffDict_),
    muC_("muC", dimless, coeffDict_),
    I0_("I0", dimless, coeffDict_),
    aQSk_("aQSk", dimless, coeffDict_),
    aInt_("aInt", dimless, coeffDict_),
    k_("k", dimensionSet(1,0,-2,0,0), coeffDict_),
    alphaDeltaMin_("alphaDeltaMin", dimless, coeffDict_),
    Rc_("Rc", dimless, coeffDict_)
{

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::frictionalStressModels::
SchneiderbauerEtAlOld::~SchneiderbauerEtAlOld()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::frictionalStressModels::
SchneiderbauerEtAlOld::frictionalPressure
(
    const phaseModel& phase,
    const dimensionedScalar& alphaMinFriction,
    const dimensionedScalar& alphaMax,
    const volScalarField& dp,
    const volScalarField& rho,
    const volSymmTensorField& D
) const
{
    const volScalarField& alpha = phase;
    volSymmTensorField S = dev(D);
        
    tmp<volScalarField> tpf
    (
        new volScalarField
        (
            IOobject
            (
                "Schneiderbauer:pf",
                phase.mesh().time().timeName(),
                phase.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            phase.mesh(),
            dimensionedScalar("pf", dimensionSet(1, -1, -2, 0, 0), 0.0)
        )
    );

    volScalarField& pf = tpf.ref();

    forAll(D, celli)
    {
        if (alpha[celli] > alphaMinFriction.value()) {
            if (alpha[celli] < alphaMax.value()) {
                pf[celli] =
                    2.0
                   *rho[celli]
                   *sqr(b_.value()*dp[celli])
                   *min(S[celli]&&S[celli],1.0e4)
                  /sqr(max(alphaMax.value() - alpha[celli],alphaDeltaMin_.value()));
            } else {
                pf[celli] =
                    aQSk_.value()
                   *k_.value()
                   *pow(max(alpha[celli] - alphaMax.value(),SMALL), 2.0/3.0)
                  /dp[celli];
            }
        }
    }

    const fvPatchList& patches = phase.mesh().boundary();
    const volVectorField& U = phase.U();

    volScalarField::Boundary& pfBf = pf.boundaryFieldRef();

    forAll(patches, patchi)
    {
        //if (!patches[patchi].coupled()) {
        const fvPatch& curPatch = patches[patchi];
        if (isA<wallFvPatch>(curPatch)) {
            pfBf[patchi] =
                pos(alpha[patchi] - alphaMinFriction.value())
               *neg(alpha[patchi] - alphaMax.value())
               *2.0
               *rho[patchi]
               *sqr(b_.value()*dp[patchi])
               *max(min(magSqr(U.boundaryField()[patchi].snGrad()),1.0e4),SMALL)
               /sqr(max(alphaMax.value() - alpha[patchi],alphaDeltaMin_.value()))
              + pos(alpha[patchi] - alphaMax.value())
               *aQSk_.value()
               *k_.value()
               *pow(max(alpha[patchi] - alphaMax.value(),SMALL), 2.0/3.0)
              /dp[patchi];
            
            Info << "max(pfW): " << max(pfBf[patchi])
                 << ", max(gradU)" << max(magSqr(U.boundaryField()[patchi].snGrad()))
                 << endl;
        }

    }

    // Correct coupled BCs
    pf.correctBoundaryConditions();

    return tpf;
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::frictionalStressModels::
SchneiderbauerEtAlOld::frictionalPressurePrime
(
    const phaseModel& phase,
    const dimensionedScalar& alphaMinFriction,
    const dimensionedScalar& alphaMax,
    const volScalarField& dp,
    const volScalarField& rho,
    const volSymmTensorField& D
) const
{
    const volScalarField& alpha = phase;
    // pPrime does not contain pInt contribution
    // pPrime is solely used, if implicitPhasePressure is true
    volSymmTensorField S = dev(D);
    volScalarField DD
    (
        min
        (
            S&&S
           ,dimensionedScalar("dmax",dimensionSet(0, 0, -2, 0, 0),1.0e4)
        )
    );
    return
         pos(alpha - alphaMinFriction)
        *neg(alpha - alphaMax)
        *4.0*rho*sqr(b_*dp)*DD
        /pow3(Foam::max(alphaMax - alpha, alphaDeltaMin_))
      +  pos(alpha- alphaMax)
        *(2.0/3.0)
        *aQSk_
        *k_
        /(pow(max(alpha - alphaMax,alphaDeltaMin_), 1.0/3.0)*dp);
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::frictionalStressModels::
SchneiderbauerEtAlOld::nu
(
    const phaseModel& phase,
    const dimensionedScalar& alphaMinFriction,
    const dimensionedScalar& alphaMax,
    const volScalarField& pf,
    const volScalarField& dp,
    const volSymmTensorField& D
) const
{
    const volScalarField& alpha = phase;

    tmp<volScalarField> tnu
    (
        new volScalarField
        (
            IOobject
            (
                "SchneiderbauerEtAlOld:nu",
                phase.mesh().time().timeName(),
                phase.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            phase.mesh(),
            dimensionedScalar("nu", dimensionSet(0, 2, -1, 0, 0), 0.0)
        )
    );

    volScalarField& nuf = tnu.ref();
    volSymmTensorField S = dev(D);

    forAll(D, celli)
    {
        if (alpha[celli] > alphaMinFriction.value())
        {
            nuf[celli] = (
                             Rc_.value()
                           + muSt_.value()
                           + (
                                muC_.value() - muSt_.value()
                             )
                            /(
                                 I0_.value()
                               / (
                                   (2.0 * sqrt(0.5*(S[celli]&&S[celli])) * dp[celli])
                                   /(sqrt(pf[celli])+SMALL) + SMALL
                                 )
                               + 1.0
                             )
                          )
                        * pf[celli]
                        / (2.0 * sqrt(0.5*(S[celli]&&S[celli])) + SMALL);
        }
    }

    const fvPatchList& patches = phase.mesh().boundary();
    const volVectorField& U = phase.U();

    volScalarField::Boundary& nufBf = nuf.boundaryFieldRef();

    forAll(patches, patchi)
    {
        // if (!patches[patchi].coupled()) {
        const fvPatch& curPatch = patches[patchi];
        if (isA<wallFvPatch>(curPatch)) {
            nufBf[patchi] = (
                               muSt_.value()
                             + (
                                   muC_.value() - muSt_.value()
                               )
                              /(
                                   I0_.value()
                                 / (
                                     (mag(U.boundaryField()[patchi].snGrad())
                                          * dp.boundaryField()[patchi])
                                     /(sqrt(pf.boundaryField()[patchi])+SMALL)
                                    + SMALL
                                   )
                                 + 1.0
                               )
                             )
                            * pf.boundaryField()[patchi]
                            / (
                                 mag(U.boundaryField()[patchi].snGrad())
                               + SMALL
                              );
        }
    }
    // Correct coupled BCs
    nuf.correctBoundaryConditions();

    return tnu;
}


bool Foam::kineticTheoryModels::frictionalStressModels::
SchneiderbauerEtAlOld::read()
{
    coeffDict_ <<= dict_.optionalSubDict(typeName + "Coeffs");

    b_.read(coeffDict_);
    muSt_.read(coeffDict_);
    muC_.read(coeffDict_);

    I0_.read(coeffDict_);
    aQSk_.read(coeffDict_);
    aInt_.read(coeffDict_);
    k_.read(coeffDict_);
    
    alphaDeltaMin_.read(coeffDict_);
    Rc_.read(coeffDict_);

    return true;
}


// ************************************************************************* //

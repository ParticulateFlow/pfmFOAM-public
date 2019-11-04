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

\*---------------------------------------------------------------------------*/

#include "SchneiderbauerEtAlFrictionalStress.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace frictionalStressModels
{
    defineTypeNameAndDebug(SchneiderbauerEtAl, 0);

    addToRunTimeSelectionTable
    (
        frictionalStressModel,
        SchneiderbauerEtAl,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::frictionalStressModels::
SchneiderbauerEtAl::SchneiderbauerEtAl
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
    aQSk_("aQSk", dimensionSet(1, 0, -2, 0, 0), coeffDict_),
    alphaDeltaMin_("alphaDeltaMin", dimless, coeffDict_)
{

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::frictionalStressModels::
SchneiderbauerEtAl::~SchneiderbauerEtAl()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::frictionalStressModels::
SchneiderbauerEtAl::frictionalPressure
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
    
    tmp<volScalarField> tpf
    (
        new volScalarField
        (
            IOobject
            (
                "SchneiderbauerEtAl:pf",
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
        if (alpha[celli] > alphaMinFriction.value())
        {
            // Equ. (32)
            pf[celli] =   2.0*rho[celli]*sqr(b_.value()*dp[celli])
                        * Foam::min(D[celli]&&D[celli],3.0e4)
                        / sqr(Foam::max(alphaMax.value() - alpha[celli], alphaDeltaMin_.value()));
        }
    }

    const fvPatchList& patches = phase.mesh().boundary();
    const volVectorField& U = phase.U();

    volScalarField::Boundary& pfBf = pf.boundaryFieldRef();

    forAll(patches, patchi)
    {
        if (!patches[patchi].coupled())
        {
            // Equ. (32)
            pfBf[patchi] =   2.0*rho.boundaryField()[patchi]*sqr(b_.value()*dp.boundaryField()[patchi])
                          * Foam::min(magSqr(U.boundaryField()[patchi].snGrad()),3.0e4)
                          / sqr(Foam::max(alphaMax.value() - alpha.boundaryField()[patchi], alphaDeltaMin_.value()));
        }
    }

    // Correct coupled BCs
    pf.correctBoundaryConditions();
    return tpf;
                
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::frictionalStressModels::
SchneiderbauerEtAl::frictionalPressurePrime
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
    
    tmp<volScalarField> tpf
    (
        new volScalarField
        (
            IOobject
            (
                "SchneiderbauerEtAl:dpf",
                phase.mesh().time().timeName(),
                phase.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            phase.mesh(),
            dimensionedScalar("dpf", dimensionSet(1, -1, -2, 0, 0), 0.0)
        )
    );

    volScalarField& pf = tpf.ref();
    
    forAll(D, celli)
    {
        if (alpha[celli] > alphaMinFriction.value())
        {
            // Equ. (32)
            pf[celli] =   4.0*rho[celli]*sqr(b_.value()*dp[celli])
                        * Foam::min(D[celli]&&D[celli],3.0e4)
                        / pow3(Foam::max(alphaMax.value() - alpha[celli], alphaDeltaMin_.value()));
        }
    }

    const fvPatchList& patches = phase.mesh().boundary();
    const volVectorField& U = phase.U();

    volScalarField::Boundary& pfBf = pf.boundaryFieldRef();

    forAll(patches, patchi)
    {
        if (!patches[patchi].coupled())
        {
            // Equ. (32)
            pfBf[patchi] =   4.0*rho.boundaryField()[patchi]*sqr(b_.value()*dp.boundaryField()[patchi])
                          * Foam::min(magSqr(U.boundaryField()[patchi].snGrad()),3.0e4)
                          / pow3(Foam::max(alphaMax.value() - alpha.boundaryField()[patchi], alphaDeltaMin_.value()));
        }
    }

    // Correct coupled BCs
    pf.correctBoundaryConditions();
    return tpf;
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::frictionalStressModels::
SchneiderbauerEtAl::nu
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
                "SchneiderbauerEtAl:nu",
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

    forAll(D, celli)
    {
        if (alpha[celli] > alphaMinFriction.value())
        {
            nuf[celli] = (
                             muSt_.value()
                           + (
                                muC_.value() - muSt_.value()
                             )
                            /(
                                 I0_.value()
                               / (
                                   (2.0 * sqrt(0.5*(D[celli]&&D[celli])) * dp[celli])
                                   /(sqrt(pf[celli])+SMALL) + SMALL
                                 )
                               + 1.0
                             )
                          )
                        * pf[celli]
                        / (2.0 * sqrt(0.5*(D[celli]&&D[celli])) + SMALL);
        }
    }

    const fvPatchList& patches = phase.mesh().boundary();
    const volVectorField& U = phase.U();

    volScalarField::Boundary& nufBf = nuf.boundaryFieldRef();

    forAll(patches, patchi)
    {
        if (!patches[patchi].coupled())
        {
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
SchneiderbauerEtAl::read()
{
    coeffDict_ <<= dict_.optionalSubDict(typeName + "Coeffs");

    b_.read(coeffDict_);
    muSt_.read(coeffDict_);
    muC_.read(coeffDict_);

    I0_.read(coeffDict_);
    aQSk_.read(coeffDict_);
    
    alphaDeltaMin_.read(coeffDict_);

    return true;
}


// ************************************************************************* //

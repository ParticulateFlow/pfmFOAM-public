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
    aQSk_("aQSk", dimless, coeffDict_),
    aInt_("aInt", dimless, coeffDict_),
    k_("k", dimensionSet(1,0,-2,0,0), coeffDict_),
    alphaDeltaMin_("alphaDeltaMin", dimless, coeffDict_),
    Rc_("Rc", dimless, coeffDict_)
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
        volScalarField DD
        (
            min
            (
                max(D&&D,dimensionedScalar("dmax",dimensionSet(0, 0, -2, 0, 0),1.0e-8))
               ,1.0e1
            )
        );
        volScalarField pInt
        (
           aInt_
           *k_
           *sqrt(sqrt(DD)*dp/sqrt(k_/(rho*dp)))
           /dp
         );

        return
             pos(alpha - alphaMinFriction)
            *neg(alpha - alphaMax)
            /(
                sqr(Foam::max(alphaMax - alpha,alphaDeltaMin_))
               /(
                    2.0
                   *rho
                   *sqr(b_*dp)
                   *DD
                )
              + 1.0
               /pInt
            )
          +
             pos(alpha - alphaMax)
            *(
                aQSk_*k_*pow(Foam::max(alpha - alphaMax, scalar(0)), 2.0/3.0)/dp
              + pInt
             );
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
        // pPrime does not contain pInt contribution
        // pPrime is solely used, if implicitPhasePressure is true
        return
             pos(alpha - alphaMinFriction)
            *4.0*rho*sqr(b_*dp)*min(D&&D,dimensionedScalar("dmax",dimensionSet(0, 0, -2, 0, 0),1.0e2))
            /pow3(Foam::max(alphaMax - alpha, alphaDeltaMin_));
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
                             Rc_.value()
                           + muSt_.value()
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
    aInt_.read(coeffDict_);
    k_.read(coeffDict_);
    
    alphaDeltaMin_.read(coeffDict_);
    Rc_.read(coeffDict_);

    return true;
}


// ************************************************************************* //

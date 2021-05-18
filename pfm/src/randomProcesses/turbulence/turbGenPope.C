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

#include "fft.H"
#include "turbGenPope.H"
#include "Kmesh.H"
#include "primitiveFields.H"
#include "Ek.H"
#include "mathematicalConstants.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::turbGenPope::turbGenPope
(
    const Kmesh& k,
    const scalar CP,
    const scalar P0,
    const scalar BETA,
    const scalar NU,
    const scalar CL,
    const scalar CETA,
    const scalar E,
    const scalar EPSILON
)
:
    K(k),
    Cp(CP),
    p0(P0),
    beta(BETA),
    nu(NU),
    cL(CL),
    cEta(CETA),
    e(E),
    epsilon(EPSILON),
    RanGen(label(0))
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::vectorField Foam::turbGenPope::U()
{
    vectorField s(K.size());
    scalarField rndPhases(K.size());

    forAll(K, i)
    {
        s[i] = RanGen.sample01<vector>();
        rndPhases[i] = RanGen.scalar01();
    }

    s = K ^ s;
    s = s/(mag(s) + 1.0e-20);

    s = Ek(Cp, p0, beta, nu, cL, cEta, e, epsilon, mag(K))*s;

    complexVectorField up
    (
        fft::reverseTransform
        (
            ComplexField(cos(constant::mathematical::twoPi*rndPhases)*s,
            sin(constant::mathematical::twoPi*rndPhases)*s),
            K.nn()
        )
    );

    return ReImSum(up);
}


// ************************************************************************* //

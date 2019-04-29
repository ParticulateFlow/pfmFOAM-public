/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "laplaceAnIsoFilter.H"
#include "addToRunTimeSelectionTable.H"
#include "calculatedFvPatchFields.H"
#include "fvm.H"
#include "fvc.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(laplaceAnIsoFilter, 0);
    addToRunTimeSelectionTable(LESfilter, laplaceAnIsoFilter, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::laplaceAnIsoFilter::laplaceAnIsoFilter(const fvMesh& mesh, scalar widthCoeff)
:
    LESfilter(mesh),
    widthCoeff_(widthCoeff),
    coeff_
    (
        IOobject
        (
            "laplaceAnIsoFilterCoeff",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("zero", dimLength*dimLength, 0),
        calculatedFvPatchScalarField::typeName
    )
{
 	volScalarField hmax
	(
    		IOobject
    		(
        		"hmax",
            		mesh.time().timeName(),
            		mesh
    		),
    		mesh,
    		dimensionedScalar("zero", dimLength, 0.0)
	);
	volScalarField deltaF
	(
    		IOobject
    		(
        		"deltaF",
            		mesh.time().timeName(),
            		mesh
    		),
    		mesh,
    		dimensionedScalar("zero", dimLength, 0.0)
	);
	volScalarField wD
	(
    		IOobject
    		(
        		"wD",
            		mesh.time().timeName(),
            		mesh
    		),
    		mesh,
    		dimensionedScalar("zero", dimLength, 0.0)
	);

    wD = wallDist(mesh).y();

	forAll(mesh.cells(),cellI)
	{
        scalar deltaMaxTmp = 0.0;
        const labelList& cFaces = mesh.cells()[cellI];
        const point& centrevector = mesh.cellCentres()[cellI];

        forAll(cFaces, cFaceI)
        {
            label faceI = cFaces[cFaceI];
            const point& facevector = mesh.faceCentres()[faceI];
            scalar tmp = mag(facevector - centrevector);
            if (tmp > deltaMaxTmp)
            {
                    deltaMaxTmp = tmp;
            }
        }
        if (neg(wD[cellI])>0.5) {
            hmax[cellI] = 2.0*deltaMaxTmp;
        } else {
            hmax[cellI] = 2.0*Foam::min(deltaMaxTmp,wD[cellI]);
        }
    }
    deltaF = hmax;
    coeff_.ref() = sqr(deltaF.ref())/widthCoeff_;
}


Foam::laplaceAnIsoFilter::laplaceAnIsoFilter(const fvMesh& mesh, const dictionary& bd)
:
    LESfilter(mesh),
    widthCoeff_(readScalar(bd.subDict(type() + "Coeffs").lookup("widthCoeff"))),
    coeff_
    (
        IOobject
        (
            "laplaceAnIsoFilterCoeff",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("zero", dimLength*dimLength, 0),
        calculatedFvPatchScalarField::typeName
    )
{

 	volScalarField hmax
	(
    		IOobject
    		(
        		"hmax",
            		mesh.time().timeName(),
            		mesh
    		),
    		mesh,
    		dimensionedScalar("zero", dimLength, 0.0)
	);
	volScalarField deltaF
	(
    		IOobject
    		(
        		"deltaF",
            		mesh.time().timeName(),
            		mesh
    		),
    		mesh,
    		dimensionedScalar("zero", dimLength, 0.0)
	);

	volScalarField wD
	(
    		IOobject
    		(
        		"wD",
            		mesh.time().timeName(),
            		mesh
    		),
    		mesh,
    		dimensionedScalar("zero", dimLength, 0.0)
	);

    wD = wallDist(mesh).y();

    forAll(mesh.cells(),cellI)
    {
        scalar deltaMaxTmp = 0.0;
        const labelList& cFaces = mesh.cells()[cellI];
        const point& centrevector = mesh.cellCentres()[cellI];
        
        forAll(cFaces, cFaceI)
        {
            label faceI = cFaces[cFaceI];
            const point& facevector = mesh.faceCentres()[faceI];
            scalar tmp = mag(facevector - centrevector);
            if (tmp > deltaMaxTmp)
            {
                deltaMaxTmp = tmp;
            }
        }
        if (neg(wD[cellI])>0.5) {
            hmax[cellI] = 2.0*deltaMaxTmp;
        } else {
            hmax[cellI] = 2.0*Foam::min(deltaMaxTmp,wD[cellI]);
        }
    }
	deltaF = hmax;
    coeff_.ref() = sqr(deltaF.ref())/widthCoeff_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::laplaceAnIsoFilter::read(const dictionary& bd)
{
    bd.subDict(type() + "Coeffs").lookup("widthCoeff") >> widthCoeff_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::laplaceAnIsoFilter::operator()
(
    const tmp<volScalarField>& unFilteredField
) const
{
    correctBoundaryConditions(unFilteredField);

    tmp<volScalarField> filteredField =
        unFilteredField() + fvc::laplacian(coeff_, unFilteredField());

    unFilteredField.clear();

    return filteredField;
}


Foam::tmp<Foam::volVectorField> Foam::laplaceAnIsoFilter::operator()
(
    const tmp<volVectorField>& unFilteredField
) const
{
    correctBoundaryConditions(unFilteredField);

    tmp<volVectorField> filteredField =
        unFilteredField() + fvc::laplacian(coeff_, unFilteredField());

    unFilteredField.clear();

    return filteredField;
}


Foam::tmp<Foam::volSymmTensorField> Foam::laplaceAnIsoFilter::operator()
(
    const tmp<volSymmTensorField>& unFilteredField
) const
{
    correctBoundaryConditions(unFilteredField);

    tmp<volSymmTensorField> filteredField =
        unFilteredField() + fvc::laplacian(coeff_, unFilteredField());

    unFilteredField.clear();

    return filteredField;
}


Foam::tmp<Foam::volTensorField> Foam::laplaceAnIsoFilter::operator()
(
    const tmp<volTensorField>& unFilteredField
) const
{
    correctBoundaryConditions(unFilteredField);

    tmp<volTensorField> filteredField =
        unFilteredField() + fvc::laplacian(coeff_, unFilteredField());

    unFilteredField.clear();

    return filteredField;
}


// ************************************************************************* //

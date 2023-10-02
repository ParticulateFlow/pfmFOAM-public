/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    computeFieldDifference


Description
    Read in two fields, take their difference and integrate over the domain

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    scalar domainVol = gSum(mesh.V());


    if (fieldType == "volScalarField")
    {
        volScalarField field0
        (
            IOobject
            (
                fieldNames[0],
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );

        volScalarField field1
        (
            IOobject
            (
                fieldNames[1],
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );

        diff.dimensions().reset(field0.dimensions());
        norm0.dimensions().reset(field0.dimensions());
        norm1.dimensions().reset(field1.dimensions());

        diff = mag(field0 - field1);
        norm0 = mag(field0);
        norm1 = mag(field1);
    }
    else if (fieldType == "volVectorField")
    {
        volVectorField field0
        (
            IOobject
            (
                fieldNames[0],
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );

        volVectorField field1
        (
            IOobject
            (
                fieldNames[1],
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );

        diff.dimensions().reset(field0.dimensions());
        norm0.dimensions().reset(field0.dimensions());
        norm1.dimensions().reset(field1.dimensions());

        diff = mag(field0 - field1);
        norm0 = mag(field0);
        norm1 = mag(field1);
    }
    else
    {
        volSymmTensorField field0
        (
            IOobject
            (
                fieldNames[0],
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );

        volSymmTensorField field1
        (
            IOobject
            (
                fieldNames[1],
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );

        diff.dimensions().reset(field0.dimensions());
        norm0.dimensions().reset(field0.dimensions());
        norm1.dimensions().reset(field1.dimensions());

        diff = mag(field0 - field1);
        norm0 = mag(field0);
        norm1 = mag(field1);
    }

    scalar integratedDiff = fvc::domainIntegrate(diff).value();
    scalar integratedMag0 = fvc::domainIntegrate(norm0).value();
    scalar integratedMag1 = fvc::domainIntegrate(norm1).value();
    
    scalar diffPerV = integratedDiff/domainVol;
    scalar diffRelative = integratedDiff / (Foam::sqrt(integratedMag0) * Foam::sqrt(integratedMag1));

    Info << "norm0 per volume = " << integratedMag0/domainVol << ", norm1 per volume = " << integratedMag1/domainVol << endl;
    Info << "difference per volume = " << diffPerV << ", relative difference = " << diffRelative << endl;

    return 0;
}

// ************************************************************************* //

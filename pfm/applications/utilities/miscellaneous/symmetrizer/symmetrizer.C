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
    symmetrizer

Description
    Read time series and symmetrize its fields (either all times or latest time)
    with respect to a given plane. Needs to be run in serial mode so that all cells
    can be found. Works only if geometry is reflection-invariant wrt the specified plane.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::addOption
    (
        "dict",
        "word",
        "specify the symmetrizerDict name"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    word dictname = "";
    if (args.optionFound("dict"))
    {
        dictname = args["dict"];
    }
    else
    {
        dictname = "symmetrizerProperties";
    }

    Info<< "dictname: " << dictname << endl;

    IOdictionary symmetrizerProperties
    (
        IOobject
        (
            dictname,
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    wordList scalarFieldNames_(symmetrizerProperties.lookup("scalarFields"));

    wordList vectorFieldNames_(symmetrizerProperties.lookup("vectorFields"));

    wordList pseudoVectorFieldNames_(symmetrizerProperties.lookup("pseudoVectorFields"));

    PtrList<volScalarField> scalarFields_;

    PtrList<volVectorField> vectorFields_;

    PtrList<volVectorField> pseudoVectorFields_;

    vector refPoint(symmetrizerProperties.lookup("refPoint"));

    vector refDirection(symmetrizerProperties.lookup("refDirection"));

    bool latestTime_(symmetrizerProperties.lookupOrDefault<bool>("latestTime",true));

    instantList timeDirs(runTime.times());

    for(int timeI = 0; timeI < timeDirs.size(); timeI++)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        if (runTime.timeName() == "constant")
        {
                continue;
        }
        if (latestTime_ && timeI != timeDirs.size()-1)
        {
                continue;
        }
        Info << "Symmetrizing fields for time = " << runTime.timeName() << endl;

        // read scalar and vector fields
        for (int i=0; i<scalarFieldNames_.size(); i++)
        {
            scalarFields_.append
            (
                new volScalarField
                (
                    IOobject
                    (
                        scalarFieldNames_[i],
                        runTime.timePath(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                     ),
                     mesh
                )
            );
        }
        for (int i=0; i<vectorFieldNames_.size(); i++)
        {
            vectorFields_.append
            (
                new volVectorField
                (
                    IOobject
                    (
                        vectorFieldNames_[i],
                        runTime.timePath(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                     ),
                     mesh
                )
            );
        }
        for (int i=0; i<pseudoVectorFieldNames_.size(); i++)
        {
            pseudoVectorFields_.append
            (
                new volVectorField
                (
                    IOobject
                    (
                        pseudoVectorFieldNames_[i],
                        runTime.timePath(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                     ),
                     mesh
                )
            );
        }

        label cellI_transformed;
        // mirror scalar fields
        for (int i=0; i<scalarFieldNames_.size(); i++)
        {
            volScalarField tmpField(scalarFields_[i]);

            forAll(tmpField, cellI)
            {
                vector position = mesh.C()[cellI];
                vector transformedPosition = 2 * ((refPoint - position) & refDirection) * refDirection / (refDirection & refDirection) + position;
                cellI_transformed = mesh.findCell(transformedPosition);
                if(cellI_transformed < 0)
                {
                    FatalError << "Couldn't find transformed cell. Aborting.\n" << abort(FatalError);
                }

                scalar value = scalarFields_[i][cellI_transformed];
                scalar transformedValue = value;
                tmpField[cellI] = 0.5 * (transformedValue + scalarFields_[i][cellI]);
            }
            scalarFields_[i] = tmpField;
            scalarFields_[i].write();
        }

        // mirror vector fields
        for (int i=0; i<vectorFieldNames_.size(); i++)
        {
            volVectorField tmpField(vectorFields_[i]);

            forAll(tmpField, cellI)
            {
                vector position = mesh.C()[cellI];
                vector transformedPosition = 2 * ((refPoint - position) & refDirection) * refDirection / (refDirection & refDirection) + position;
                cellI_transformed = mesh.findCell(transformedPosition);
                if(cellI_transformed < 0)
                {
                    FatalError << "Couldn't find transformed cell. Aborting.\n" << abort(FatalError);
                }

                vector value = vectorFields_[i][cellI_transformed];
                vector transformedValue = -2 * (value & refDirection) * refDirection / (refDirection & refDirection) + value;
                tmpField[cellI] = 0.5 * (transformedValue + vectorFields_[i][cellI]);
            }
            vectorFields_[i] = tmpField;
            vectorFields_[i].write();
        }

        // mirror pseudovector fields
        for (int i=0; i<pseudoVectorFieldNames_.size(); i++)
        {
            volVectorField tmpField(pseudoVectorFields_[i]);

            forAll(tmpField, cellI)
            {
                vector position = mesh.C()[cellI];
                vector transformedPosition = 2 * ((refPoint - position) & refDirection) * refDirection / (refDirection & refDirection) + position;
                cellI_transformed = mesh.findCell(transformedPosition);
                if(cellI_transformed < 0)
                {
                    FatalError << "Couldn't find transformed cell. Aborting.\n" << abort(FatalError);
                }

                vector value = pseudoVectorFields_[i][cellI_transformed];
                vector transformedValue = value;
                tmpField[cellI] = 0.5 * (transformedValue + pseudoVectorFields_[i][cellI]);
            }
            pseudoVectorFields_[i] = tmpField;
            pseudoVectorFields_[i].write();
        }


        scalarFields_.clear();
        vectorFields_.clear();
        pseudoVectorFields_.clear();
    }


    Info << "\nEnd" << endl;

    return 0;
}

// ************************************************************************* //

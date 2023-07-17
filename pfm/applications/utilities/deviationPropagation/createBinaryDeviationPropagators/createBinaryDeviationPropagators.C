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

Application
    createBinaryDeviationPropagators

Description
    Reads all deviation propagators of a time series stored in ASCII format and
    saves them in binary format.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "IFstream.H"
#include "OFstream.H"
#include "tensorList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    IOdictionary dataBaseProperties
    (
        IOobject
        (
            "dataBaseProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    fileNameList dataBases(dataBaseProperties.lookupOrDefault<fileNameList>("dataBases", fileNameList(1,"dataBase")));

 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    if (!Pstream::master())
    {
        FatalError <<"current implementation for serial runs only\n" << abort(FatalError);
    }

    forAll(dataBases,i)
    {
        fileName dbName = dataBases[i];
        Info << "\nReading deviation propagators of database " << dbName << endl;

        Foam::Time dbTime(dbName, "", "../system", "../constant", false);
        instantList timeDirs(dbTime.times());
        if (timeDirs.size() == 0)
        {
            FatalError <<"database " << dbName << " does not exist or is empty\n" << abort(FatalError);
        }

        for (instantList::iterator it=timeDirs.begin(); it != timeDirs.end(); ++it)
        {
            // set time
            dbTime.setTime(*it, it->value());

            // skip constant
            if (dbTime.timeName() == "constant")
            {
                continue;
            }

            Info << "Reading at t = " << dbTime.timeName() << endl;

            List<tensorList> K_internal;
            IFstream IS_internal(dbName+"/"+dbTime.timeName()+"/K_uu_internal.txt");
            IS_internal >> K_internal;
            OFstream OS_internal(dbName+"/"+dbTime.timeName()+"/K_uu_internal", IOstream::BINARY);
            OS_internal << K_internal;

            List<tensorList> K_boundary;
            IFstream IS_boundary(dbName+"/"+dbTime.timeName()+"/K_uu_boundary.txt");
            IS_boundary >> K_boundary;
            OFstream OS_boundary(dbName+"/"+dbTime.timeName()+"/K_uu_boundary", IOstream::BINARY);
            OS_boundary << K_boundary;

            tensorList K_integrated;
            IFstream IS_integrated(dbName+"/"+dbTime.timeName()+"/K_uu_integrated.txt");
            IS_integrated >> K_integrated;
            OFstream OS_integrated(dbName+"/"+dbTime.timeName()+"/K_uu_integrated", IOstream::BINARY);
            OS_integrated << K_integrated;
        }
        Info << "Conversion of deviation propagators of database " << dbName <<" done\n" << endl;
    }



    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //


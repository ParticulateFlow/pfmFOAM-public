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

#include "sampleCellSets.H"
#include "fieldTypes.H"
#include "addToRunTimeSelectionTable.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(sampleCellSets, 0);
    addToRunTimeSelectionTable(functionObject, sampleCellSets, dictionary);
}
}



// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::sampleCellSets::writeFileHeader(const label i)
{
 /*   OFstream& file = this->file();

    writeHeader(file, "Field minima and maxima");
    writeCommented(file, "Time");

    if (location_)
    {
        writeTabbed(file, "field");

        writeTabbed(file, "min");
        writeTabbed(file, "location(min)");

        if (Pstream::parRun())
        {
            writeTabbed(file, "processor");
        }

        writeTabbed(file, "max");
        writeTabbed(file, "location(max)");

        if (Pstream::parRun())
        {
            writeTabbed(file, "processor");
        }
    }
    else
    {
        forAll(fieldSet_, fieldi)
        {
            writeTabbed(file, "min(" + fieldSet_[fieldi] + ')');
            writeTabbed(file, "max(" + fieldSet_[fieldi] + ')');
        }
    }

    file<< endl;
*/
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::sampleCellSets::sampleCellSets
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObject(name),
    distDataMeshName_(dict.lookupOrDefault<word>("distDataMesh","region0")),
    distDataMesh_(refCast<const fvMesh>(runTime.lookupObject<objectRegistry>(distDataMeshName_))),
    cellSetNames_(),
    distDataScalarFieldNames_(),
    distDataVectorFieldNames_(),
    location_(true)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::sampleCellSets::~sampleCellSets()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::sampleCellSets::read(const dictionary& dict)
{
    location_ = dict.lookupOrDefault<Switch>("location", true);


    dict.lookup("distDataScalarFields") >> distDataScalarFieldNames_;

    dict.lookup("distDataVectorFields") >> distDataVectorFieldNames_;
    
    dict.lookup("cellSets") >> cellSetNames_;

    numCellSets_ = cellSetNames_.size();

    for (label i=0; i<numCellSets_; i++)
    {
        cellSets_.append(new cellSet(distDataMesh_,cellSetNames_[i]));
    }

    return true;
}


bool Foam::functionObjects::sampleCellSets::execute()
{
    label nScalarFields = distDataScalarFieldNames_.size();
    label nVectorFields = distDataVectorFieldNames_.size();

    // if time directories do not exist yet, create them
    mkDir(distDataMesh_.time().timeName());

    // if time to write fields values in cell set to textfiles, do it
    for (label cellset = 0; cellset < numCellSets_; cellset++)
    {
        fileName filename = distDataMesh_.time().timeName()+"/data_"+cellSetNames_[cellset];
        OFstream sampleFile(filename);
        labelList cellsInSet = cellSets_[cellset].toc();
        scalarRectangularMatrix dataInCellSet(cellsInSet.size(),nScalarFields+3*nVectorFields);

        sampleFile << "# c_x\tc_y\tc_z\t";
        // first sample scalar fields
        for (label fieldI = 0; fieldI < nScalarFields; fieldI++)
        {
            sampleFile << distDataScalarFieldNames_[fieldI] << "\t";
            const volScalarField& field = distDataMesh_.lookupObject<volScalarField>(distDataScalarFieldNames_[fieldI]);
            forAll(cellsInSet, cellI)
            {
                label cellID = cellsInSet[cellI];
                dataInCellSet(cellI,fieldI) = field[cellID];
            }
        }

        // first sample vector fields
        for (label fieldI = 0; fieldI < nVectorFields; fieldI++)
        {
            sampleFile << distDataVectorFieldNames_[fieldI] << "_x\t";
            sampleFile << distDataVectorFieldNames_[fieldI] << "_y\t";
            sampleFile << distDataVectorFieldNames_[fieldI] << "_z\t";
            const volVectorField& field = distDataMesh_.lookupObject<volVectorField>(distDataVectorFieldNames_[fieldI]);
            forAll(cellsInSet, cellI)
            {
                label cellID = cellsInSet[cellI];
                dataInCellSet(cellI,nScalarFields+3*fieldI) = field[cellID].x();
                dataInCellSet(cellI,nScalarFields+3*fieldI+1) = field[cellID].y();
                dataInCellSet(cellI,nScalarFields+3*fieldI+2) = field[cellID].z();
            }

        }

        sampleFile << endl << "#";
        for (label i=0; i<nScalarFields + 3*nVectorFields+3; i++)
        {
            sampleFile << "---------";
        }
        sampleFile << endl;

        // now print data matrix
        for (label row=0; row<dataInCellSet.m(); row++)
        {
            // first print cell coordinates, then field values
            label cellID = cellsInSet[row];
            scalar cx = distDataMesh_.C()[cellID].x();
            scalar cy = distDataMesh_.C()[cellID].y();
            scalar cz = distDataMesh_.C()[cellID].z();
            sampleFile <<  cx << "\t" << cy << "\t" << cz << "\t";
            for (label col=0; col<dataInCellSet.n(); col++)
            {
                sampleFile << dataInCellSet(row,col) << "\t";
            }
            sampleFile << endl;
        }

    }
    return true;
}


bool Foam::functionObjects::sampleCellSets::write()
{

    return true;
}



// ************************************************************************* //

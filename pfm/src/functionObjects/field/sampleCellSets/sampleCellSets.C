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
    distDataMeshName_(dict.lookup("distDataMesh")),
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
    // if time to write fields values in cell set to textfiles, do it
    for (label cellset = 0; cellset < numCellSets_; cellset++)
    {
        // first sample scalar fields
        for (label fieldI = 0; fieldI < distDataScalarFieldNames_.size(); fieldI++)
        {
            fileName filename = distDataMesh_.time().timeName()+"/data_"+cellSetNames_[cellset]+"_"+distDataScalarFieldNames_[fieldI];
            OFstream sampleFile(filename);

        }

        // first sample vector fields
        for (label fieldI = 0; fieldI < distDataVectorFieldNames_.size(); fieldI++)
        {

        }
    }
    return true;
}


bool Foam::functionObjects::sampleCellSets::write()
{

    return true;
}



// ************************************************************************* //

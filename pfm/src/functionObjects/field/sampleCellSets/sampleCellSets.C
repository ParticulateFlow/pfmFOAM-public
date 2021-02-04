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
    dx_(0.0),
    dy_(0.0),
    dz_(0.0),
    distDataScalarFieldNames_(),
    distDataVectorFieldNames_(),
    timeListFile_("sampleTimes")
{
    read(dict);

    // get cell sets extensions and write them to file for later checks
    OFstream cellSetsFile("cellSetsExtensions");
    cellSetsFile << "# cellSetIndex\txmin\txmax\tymin\tymax\tzmin\tzmax\tnx\tny\tnz" << endl;
    for (label cellset = 0; cellset < numCellSets_; cellset++)
    {
        scalar xmin = 1e6;
        scalar xmax = -1e6;
        scalar ymin = 1e6;
        scalar ymax = -1e6;
        scalar zmin = 1e6;
        scalar zmax = -1e6;
        labelList cellsInSet = cellSets_[cellset].toc();
        // get min and max cell coordinates
        forAll(cellsInSet, cellI)
        {
            label cellID = cellsInSet[cellI];
            scalar cx = distDataMesh_.C()[cellID].x();
            scalar cy = distDataMesh_.C()[cellID].y();
            scalar cz = distDataMesh_.C()[cellID].z();
            if (cx > xmax) xmax = cx;
            else if (cx < xmin) xmin = cx;
            if (cy > ymax) ymax = cy;
            else if (cy < ymin) ymin = cy;
            if (cz > zmax) zmax = cz;
            else if (cz < zmin) zmin = cz;
        }

        // get next-to-max coordinates
        scalar xmax2 = -1e6;
        scalar ymax2 = -1e6;
        scalar zmax2 = -1e6;
        forAll(cellsInSet, cellI)
        {
            label cellID = cellsInSet[cellI];
            scalar cx = distDataMesh_.C()[cellID].x();
            scalar cy = distDataMesh_.C()[cellID].y();
            scalar cz = distDataMesh_.C()[cellID].z();
            if (cx > xmax2 && cx < xmax-1e-10) xmax2 = cx;
            if (cy > ymax2 && cy < ymax-1e-10) ymax2 = cy;
            if (cz > zmax2 && cz < zmax-1e-10) zmax2 = cz;
        }

        label nx = 1;
        label ny = 1;
        label nz = 1;
        // check if box is 2D
        if (xmax2 > -9e5)
        {
            dx_ = xmax - xmax2;
            nx = round((xmax - xmin) / dx_) + 1;
        }
        if (ymax2 > -9e5)
        {
            dy_ = ymax - ymax2;
            ny = round((ymax - ymin) / dy_) + 1;
        }
        if (zmax2 > -9e5)
        {
            dz_ = zmax - zmax2;
            nz = round((zmax - zmin) / dz_) + 1;
        }

        boxExtentsN_[cellset].append(nx);
        boxExtentsN_[cellset].append(ny);
        boxExtentsN_[cellset].append(nz);

        cellSetsFile << cellset << "\t\t" << xmin << "\t" << xmax << "\t" << ymin << "\t" << ymax <<
                        "\t" << zmin << "\t" << zmax << "\t" << nx << "\t" << ny << "\t" << nz << endl;
        boxExtents_[cellset].append(xmin);
        boxExtents_[cellset].append(xmax);
        boxExtents_[cellset].append(ymin);
        boxExtents_[cellset].append(ymax);
        boxExtents_[cellset].append(zmin);
        boxExtents_[cellset].append(zmax);
    }

    // consistency check for box-shaped sets with equal cell and box sizes
    for (label cellset = 0; cellset < numCellSets_; cellset++)
    {
        label nx = boxExtentsN_[cellset][0];
        label ny = boxExtentsN_[cellset][1];
        label nz = boxExtentsN_[cellset][2];
        labelList cellsInSet = cellSets_[cellset].toc();
        if (nx * ny * nz != cellsInSet.size())
        {
            FatalError << "consistency check for box-shape of cell sets failed for set " << cellSetNames_[cellset] << 
                ": nx = " << nx << ", ny = " << ny << ", nz = " << nz << ", number of cells = " << cellsInSet.size() << abort(FatalError);
        }
        if (cellset < numCellSets_-1)
        {
            if (nx != boxExtentsN_[cellset+1][0] || ny != boxExtentsN_[cellset+1][1] || nz != boxExtentsN_[cellset+1][2])
            {
                FatalError << "consistency check for congruence of cell sets failed for sets " << cellSetNames_[cellset] << 
                    " and " << cellSetNames_[cellset+1] << abort(FatalError);
            }
        }
    }

    // write out header file for order of field output
    fileName filenameHeader = "distData_header";
    OFstream sampleFileHeader(filenameHeader);
    labelList cellsInSet = cellSets_[0].toc();
    label numCells = cellsInSet.size();
    label nScalarFields = distDataScalarFieldNames_.size();
    label nVectorFields = distDataVectorFieldNames_.size();

    sampleFileHeader << "# ";
    for (label fieldI = 0; fieldI < nScalarFields; fieldI++)
    {
        word fieldname = distDataScalarFieldNames_[fieldI];
        for (label cell=0; cell<numCells; cell++)
        {
            sampleFileHeader << fieldname << "_" << cell; 
            if (distDataVectorFieldNames_.size() > 0 || cell < numCells-1)
            {
                sampleFileHeader << ", ";
            }
        }
    }
    for (label fieldI = 0; fieldI < nVectorFields; fieldI++)
    {
        word fieldname = distDataVectorFieldNames_[fieldI];
        for (label cell=0; cell<numCells; cell++)
        {
            sampleFileHeader << fieldname << "X_" << cell << ", "; 
        }

        fieldname = distDataVectorFieldNames_[fieldI];
        for (label cell=0; cell<numCells; cell++)
        {
            sampleFileHeader << fieldname << "Y_" << cell << ", "; 
        }

        fieldname = distDataVectorFieldNames_[fieldI];
        for (label cell=0; cell<numCells; cell++)
        {
            sampleFileHeader << fieldname << "Z_" << cell;
            if (fieldI < nVectorFields - 1 || cell < numCells - 1)
            {
                sampleFileHeader << ", ";
            }
        }
    }  
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::sampleCellSets::~sampleCellSets()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::sampleCellSets::read(const dictionary& dict)
{
    dict.lookup("distDataScalarFields") >> distDataScalarFieldNames_;

    dict.lookup("distDataVectorFields") >> distDataVectorFieldNames_;
    
    dict.lookup("cellSets") >> cellSetNames_;

    numCellSets_ = cellSetNames_.size();

    boxExtents_.setSize(numCellSets_);

    boxExtentsN_.setSize(numCellSets_);

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

    timeListFile_ << distDataMesh_.time().value() << endl;

    // if time to write fields values in cell set to textfiles, do it
    for (label cellset = 0; cellset < numCellSets_; cellset++)
    {
        fileName filenameMat = distDataMesh_.time().timeName()+"/distData_cellSet_"+Foam::name(cellset)+"_matrix";
        OFstream sampleFileMat(filenameMat);
        fileName filenameCSV = distDataMesh_.time().timeName()+"/distData_cellSet_"+Foam::name(cellset)+".csv";
        OFstream sampleFileCSV(filenameCSV);
        labelList cellsInSet = cellSets_[cellset].toc();
        scalarRectangularMatrix dataInCellSet(cellsInSet.size(),nScalarFields+3*nVectorFields);

        sampleFileMat << "# cX\tcY\tcZ\t";
        // first sample scalar fields
        for (label fieldI = 0; fieldI < nScalarFields; fieldI++)
        {
            sampleFileMat << distDataScalarFieldNames_[fieldI] << "\t";
            const volScalarField& field = distDataMesh_.lookupObject<volScalarField>(distDataScalarFieldNames_[fieldI]);
            vector pos;
            label cellID;
            label counter = 0;
            // make sure output for all sets is ordered the same way (first runs z, then y, then x)
            for (label ix = 0; ix < boxExtentsN_[cellset][0]; ix++)
            {
                for (label iy = 0; iy < boxExtentsN_[cellset][1]; iy++)
                {
                    for (label iz = 0; iz < boxExtentsN_[cellset][2]; iz++)
                    {
                        pos.x() = boxExtents_[cellset][0] + ix*dx_;
                        pos.y() = boxExtents_[cellset][2] + iy*dy_;
                        pos.z() = boxExtents_[cellset][4] + iz*dz_;
                        cellID = distDataMesh_.findCell(pos);
                        dataInCellSet(counter,fieldI) = field[cellID];
                        counter++;
                    }
                }
            }
        }

        // first sample vector fields
        for (label fieldI = 0; fieldI < nVectorFields; fieldI++)
        {
            sampleFileMat << distDataVectorFieldNames_[fieldI] << "X\t";
            sampleFileMat << distDataVectorFieldNames_[fieldI] << "Y\t";
            sampleFileMat << distDataVectorFieldNames_[fieldI] << "Z\t";
            const volVectorField& field = distDataMesh_.lookupObject<volVectorField>(distDataVectorFieldNames_[fieldI]);
            vector pos;
            label cellID;
            label counter = 0;
            // make sure output for all sets is ordered the same way (first runs z, then y, then x)
            for (label ix = 0; ix < boxExtentsN_[cellset][0]; ix++)
            {
                for (label iy = 0; iy < boxExtentsN_[cellset][1]; iy++)
                {
                    for (label iz = 0; iz < boxExtentsN_[cellset][2]; iz++)
                    {
                        pos.x() = boxExtents_[cellset][0] + ix*dx_;
                        pos.y() = boxExtents_[cellset][2] + iy*dy_;
                        pos.z() = boxExtents_[cellset][4] + iz*dz_;
                        cellID = distDataMesh_.findCell(pos);
                        dataInCellSet(counter,nScalarFields+3*fieldI) = field[cellID].x();
                        dataInCellSet(counter,nScalarFields+3*fieldI+1) = field[cellID].y();
                        dataInCellSet(counter,nScalarFields+3*fieldI+2) = field[cellID].z();
                        counter++;
                    }
                }
            }
        }

        sampleFileMat << endl << "#";
        for (label i=0; i<nScalarFields + 3*nVectorFields+3; i++)
        {
            sampleFileMat << "---------";
        }
        sampleFileMat << endl;

        // print data in matrix format
        for (label row=0; row<dataInCellSet.m(); row++)
        {
            // first print cell coordinates, then field values
            label cellID = cellsInSet[row];
            scalar cx = distDataMesh_.C()[cellID].x();
            scalar cy = distDataMesh_.C()[cellID].y();
            scalar cz = distDataMesh_.C()[cellID].z();
            sampleFileMat <<  cx << "\t" << cy << "\t" << cz << "\t";
            for (label col=0; col<dataInCellSet.n(); col++)
            {
                sampleFileMat << dataInCellSet(row,col) << "\t";
            }
            sampleFileMat << endl;
        }

        // print data as csv file
        for (label col=0; col<dataInCellSet.n(); col++)
        {
            for (label row=0; row<dataInCellSet.m(); row++)
            {
                sampleFileCSV << dataInCellSet(row,col); 
                if (col < dataInCellSet.n()-1 || row < dataInCellSet.m()-1)
                {
                    sampleFileCSV << ", ";
                }
            }
        }
    }
    return true;
}


bool Foam::functionObjects::sampleCellSets::write()
{
    return true;
}

label Foam::functionObjects::sampleCellSets::round(scalar x)
{
    int ix = static_cast <int> (x);
    if (x - ix > 0.5) ix++;
    return ix;
}

// ************************************************************************* //

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

Class
    Foam::functionObjects::sampleCellSets

Description
    Calculates the value and location of scalar minimum and maximum for a list
    of user-specified fields.

    For variables with a rank greater than zero, either the min/max of a
    component value or the magnitude is reported.  When operating in parallel,
    the processor owning the value is also given.

    Example of function object specification:
    \verbatim
    sampleCellSets1
    {
        type        sampleCellSets;
        libs        ("libfieldFunctionObjects.so");
        ...
        write       yes;
        log         yes;
        location    yes;
        mode        magnitude;
        fields
        (
            U
            p
        );
    }
    \endverbatim

Usage
    \table
        Property     | Description             | Required    | Default value
        type         | type name: sampleCellSets  | yes         |
        write        | write min/max data to file |  no      | yes
        log          | write min/max data to standard output | no | no
        location     | write location of the min/max value | no | yes
        mode         | calculation mode: magnitude or component | no | magnitude
    \endtable

    Output data is written to the file \<timeDir\>/sampleCellSets.dat

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::logFiles

SourceFiles
    sampleCellSets.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_sampleCellSets_H
#define functionObjects_sampleCellSets_H

#include "vector.H"
#include "cellSet.H"
#include "fvCFD.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class sampleCellSets Declaration
\*---------------------------------------------------------------------------*/

class sampleCellSets:
    public functionObject
{
public:



protected:

    // Protected data

    word distDataMeshName_;

    const fvMesh& distDataMesh_;

    wordList cellSetNames_;

    label numCellSets_;

    PtrList<cellSet> cellSets_;

    List<scalarList> boxExtents_;

    List<labelList> boxExtentsN_;

    scalar dx_;

    scalar dy_;

    scalar dz_;

    wordList distDataScalarFieldNames_;

    wordList distDataVectorFieldNames_;

    OFstream timeListFile_;

    label round(scalar);

    // Protected Member Functions

  
        //- Output file header information
        virtual void writeFileHeader(const label i);


public:

    //- Runtime type information
    TypeName("sampleCellSets");


    // Constructors

        //- Construct from Time and dictionary
        sampleCellSets
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~sampleCellSets();


    // Member Functions

        //- Read the field min/max data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the sampleCellSets
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "sampleCellSetsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

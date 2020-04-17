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

Class
    Foam::functionObjects::mapAveField

Group
    grpFieldFunctionObjects

Description
    This function object calculates average quantities for a user-specified
    selection of volumetric fields using fieldAverage class and map the averaged
    fields on to a coarser/finer mesh.

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_mapAveField_H
#define functionObjects_mapAveField_H

#include "fvMeshFunctionObject.H"
#include "fieldAverage.H"
#include "meshToMesh.H"
#include "fvMesh.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;
class myClass;

namespace functionObjects
{

// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                        Class mapAveField Declaration
\*---------------------------------------------------------------------------*/

class mapAveField
 :
   public fieldAverage,
   public fvMesh
{

 protected:

    bool mappingInitialised_;

    const objectRegistry& obrTgt_;
    const fvMesh& rMesh_;

    labelPairList facePairs_;

    meshToMesh* cfdTorcfdInterpPtr;

    void findPais();

    void addMappedfields();

    template<class Type,class Type2>
    void addMappedFieldType(const label fieldi);

    template<class Type>
    void addMappedField(const label fieldi);

    template<class Type1, class Type2>
    void addMeanFieldSqrType(const label fieldi);

    template<class Type1, class Type2>
    void addMeanFieldSqr(const label fieldi);

    template<class Type1, class Type2>
    void addMappedPrimeFieldType(const label fieldi);

    template<class Type1, class Type2>
    void addMappedPrimeField(const label fieldi);

    template<class Type1, class Type2>
    void calcMeanFieldSqrType(const label fieldi);

    template<class Type1, class Type2>
    void calcMeanFieldSqr();

    virtual void mapAveragedField();

    template<class Type,class Type2>
    void mapMeanFieldType(const label fieldi);

    template<class Type, class Type2>
    void mapMeanSFieldType(const label fieldi);

    template<class Type>
    void mapMeanFields();

    template<class Type1, class Type2, class Type>
    void mapPrime2MeanFieldType(const label fieldi);

    template<class Type1, class Type2>
    void mapPrime2MeanField();

    void calculateAverageSqr();

    virtual void writeMappedAverages();

    void writeMappedFields() const;

    //- Write fields
    template<class Type>
    void writeMappedFieldType(const word& fieldName) const;

    //- Write fields
    template<class Type>
    void writeMappedField() const;

 public:

    //- Runtime type information
    TypeName("mapAveField");


    // Constructors

        //- Construct from Time and dictionary
        mapAveField
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );

    //- Destructor
    virtual ~mapAveField();

    //-member functions
      //virtual bool calc();
      virtual bool read(const dictionary&);
      virtual bool execute();
      virtual bool write();
};

}
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "mapAveFieldTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

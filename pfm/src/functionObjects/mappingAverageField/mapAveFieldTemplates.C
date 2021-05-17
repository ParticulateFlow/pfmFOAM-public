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

#include "fieldAverageItem.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "OFstream.H"
#include "mapAveField.H"
#include "meshToMesh.H"
#include "fvMesh.H"
#include "fvCFD.H"
#include "autoPtr.H"
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//========================= adding fields ===================================//
template<class Type,class Type2>
void Foam::functionObjects::mapAveField::addMappedFieldType(const label fieldi)
{
    const word& mappedFieldName = "mapped"+faItems_[fieldi].meanFieldName();
    const word& fieldName = faItems_[fieldi].meanFieldName();

    const Type& baseField = obr_.lookupObject<Type>(fieldName);

    Log << "    Reading/adding mapped field " << mappedFieldName << "  with the type of " << Type::typeName << endl;

        // Store on registry
        obrTgt_.store
        (
            new Type
            (
                IOobject
                (
                    mappedFieldName,
                    obrTgt_.time().timeName(obrTgt_.time().startTime().value()),
                    obrTgt_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                rMesh_,
                dimensioned<Type2>(mappedFieldName,baseField.dimensions(),pTraits<Type2>::zero)
            )
        );
}

template<class Type>
void Foam::functionObjects::mapAveField::addMappedField(const label fieldi)
{
    if (faItems_[fieldi].mean())
    {
        typedef GeometricField<Type, fvPatchField, volMesh>VolFieldType;
        typedef GeometricField<Type, fvsPatchField, surfaceMesh>SurfaceFieldType;
        const word& fieldName = faItems_[fieldi].meanFieldName();

        if (obr_.foundObject<VolFieldType>(fieldName))
        {
            addMappedFieldType<VolFieldType,Type>(fieldi);
        }
        else if (obr_.foundObject<SurfaceFieldType>(fieldName))
        {
            addMappedFieldType<SurfaceFieldType,Type>(fieldi);
        }
    }
}

template<class Type1, class Type2>
void Foam::functionObjects::mapAveField::addMeanFieldSqrType(const label fieldi)
{
    const word& fieldName = faItems_[fieldi].fieldName();
    const Type1& baseField = obr_.lookupObject<Type1>(fieldName);
    const word& fieldSqrMean = faItems_[fieldi].fieldName()+"SqrMean";

    obr_.store
    (
        new Type2
        (
            IOobject
            (
                fieldSqrMean,
                obr_.time().timeName(obr_.time().startTime().value()),
                obr_,
                restartOnOutput_
              ? IOobject::NO_READ
              : IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            magSqr(baseField)
        )
    );

}

template<class Type1, class Type2>
void Foam::functionObjects::mapAveField::addMeanFieldSqr(const label fieldi)
{
    typedef GeometricField<Type1, fvPatchField, volMesh> VolFieldType1;
    typedef GeometricField<Type2, fvPatchField, volMesh> VolFieldType2;

    if(faItems_[fieldi].prime2Mean())
    {
        const word& fieldName = faItems_[fieldi].fieldName();
        if (obr_.foundObject<VolFieldType1>(fieldName))
        {
            addMeanFieldSqrType<VolFieldType1, VolFieldType2>(fieldi);
        }

    }
}

template<class Type1, class Type2>
void Foam::functionObjects::mapAveField::addMappedPrimeFieldType(const label fieldi)
{
    const word& mappedFieldName = "mapped"+faItems_[fieldi].meanFieldName();
    const Type1& mappedMeanField = obrTgt_.lookupObject<Type1>(mappedFieldName);
    const word& mappedPrimeField = "mapped"+faItems_[fieldi].prime2MeanFieldName();

    obrTgt_.store
    (
        new Type2
        (
            IOobject
            (
                mappedPrimeField,
                obrTgt_.time().timeName(obrTgt_.time().startTime().value()),
                obrTgt_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            //sqr(mappedMeanField)
            magSqr(mappedMeanField)
        )
    );
}

template<class Type1, class Type2>
void Foam::functionObjects::mapAveField::addMappedPrimeField(const label fieldi)
{
    typedef GeometricField<Type1, fvPatchField, volMesh> VolFieldType1;
    typedef GeometricField<Type2, fvPatchField, volMesh> VolFieldType2;

    if(faItems_[fieldi].prime2Mean())
    {
        const word& fieldName = faItems_[fieldi].fieldName();
        if (obr_.foundObject<VolFieldType1>(fieldName))
        {
            addMappedPrimeFieldType<VolFieldType1, VolFieldType2>(fieldi);
        }
    }
}

//========================= square fields ===================================//
template<class Type1, class Type2>
void Foam::functionObjects::mapAveField::calcMeanFieldSqrType(const label fieldi)
{
    const word& fieldName = faItems_[fieldi].fieldName();
    const Type1& baseField = obr_.lookupObject<Type1>(fieldName);

    Type2& fieldSqrMean = const_cast<Type2&>
      (
            obr_.lookupObject<Type2>(faItems_[fieldi].fieldName()+"SqrMean")
      );

    scalar dt = obr_.time().deltaTValue();
    scalar Dt = totalTime_[fieldi];

    Dt = Dt - dt;

    if (faItems_[fieldi].iterBase())
    {
        dt = 1.0;
        Dt = scalar(totalIter_[fieldi]);
        Dt = Dt - 1;
    }

    scalar alpha = (Dt - dt)/Dt;
    scalar beta = dt/Dt;

    if (faItems_[fieldi].window() > 0)
    {
        const scalar w = faItems_[fieldi].window();

        if (Dt - dt >= w)
        {
            alpha = (w - dt)/w;
            beta = dt/w;
        }
    }

    //fieldSqrMean = alpha*fieldSqrMean + beta*sqr(baseField);
    fieldSqrMean = alpha*fieldSqrMean + beta*magSqr(baseField);

}

template<class Type1, class Type2>
void Foam::functionObjects::mapAveField::calcMeanFieldSqr()
{
    typedef GeometricField<Type1, fvPatchField, volMesh> VolFieldType1;
    typedef GeometricField<Type2, fvPatchField, volMesh> VolFieldType2;

     forAll(faItems_, fieldi)
     {
         if (faItems_[fieldi].prime2Mean())
         {
             const word& fieldName = faItems_[fieldi].fieldName();
             if (obr_.foundObject<VolFieldType1>(fieldName))
             {
                 calcMeanFieldSqrType<VolFieldType1, VolFieldType2>(fieldi);
             }
         }
     }
}

//========================= mapping fields ===================================//

template<class Type, class Type2>
void Foam::functionObjects::mapAveField::mapMeanFieldType
(
    const label fieldi
)
{

    const word& fieldName = faItems_[fieldi].meanFieldName();

    const Type& baseField = obr_.lookupObject<Type>(fieldName);

     Type& mappedField = const_cast<Type&>
      (
            obrTgt_.lookupObject<Type>("mapped"+faItems_[fieldi].meanFieldName())
      );

     cfdTorcfdInterpPtr->mapSrcToTgt(baseField,plusEqOp<Type2>(),mappedField);

    Log << " End of the mapping mean field" << endl;

}

template<class Type, class Type2>
void Foam::functionObjects::mapAveField::mapMeanSFieldType
(
    const label fieldi
)
{

    const word& fieldName = faItems_[fieldi].meanFieldName();

    const Type& baseField = obr_.lookupObject<Type>(fieldName);

    const fvMesh& mesh_ = refCast<const fvMesh>(obr_);

    Type& mappedField = const_cast<Type&>
      (
            obrTgt_.lookupObject<Type>("mapped"+faItems_[fieldi].meanFieldName())
      );

         forAll(facePairs_,i)
         {

             // check if face is an internal face (yes: returns true, no: false)
             bool internalFace = rMesh_.isInternalFace(facePairs_[i].second());

             // check if face is located on boundary patch (yes: returns boundary ID, no: -1)
             label srcPatchID = mesh_.boundaryMesh().whichPatch(facePairs_[i].first());
             label tgtPatchID = rMesh_.boundaryMesh().whichPatch(facePairs_[i].second());

             if(internalFace)
             {
                 mappedField[facePairs_[i].second()] += baseField[facePairs_[i].first()];

                 Log << " source_face_internal " << facePairs_[i].first() << ", value = " << baseField[facePairs_[i].first()] << endl;
                 Log << " target_face_internal " << facePairs_[i].second() << ", value = " << mappedField[facePairs_[i].second()] << endl;
             }
             else if(tgtPatchID > -1)
             {
                 // access field values depending on face being on boundary or internal
                 label srcFaceID = mesh_.boundaryMesh()[srcPatchID].whichFace(facePairs_[i].first());
                 label tgtFaceID = rMesh_.boundaryMesh()[tgtPatchID].whichFace(facePairs_[i].second());

                 mappedField.boundaryFieldRef()[tgtPatchID][tgtFaceID] += baseField.boundaryField()[srcPatchID][srcFaceID];

                 Log << " source_face_boundary " << facePairs_[i].first() << ", value = " << baseField.boundaryField()[srcPatchID][srcFaceID] << endl;
                 Log << " target_face_boundary " << facePairs_[i].second() << ", value = " << mappedField.boundaryField()[tgtPatchID][tgtFaceID] << endl;
             }

         }

    Log << " End of the mapping mean surface field" << endl;

}

template<class Type>
void Foam::functionObjects::mapAveField::mapMeanFields()
{

    forAll(faItems_, i)
    {
        if (faItems_[i].mean())
        {
            typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
            typedef GeometricField<Type, fvsPatchField, surfaceMesh>SurfaceFieldType;
            const word& fieldName = faItems_[i].meanFieldName();

            if (obr_.foundObject<VolFieldType>(fieldName))
            {
                Log << " Mapping the field of : " << faItems_[i].meanFieldName() << "  with the type of " << VolFieldType::typeName << nl;
                mapMeanFieldType<VolFieldType,Type>(i);
            }
            else if (obr_.foundObject<SurfaceFieldType>(fieldName))
            {
                Log << " Mapping the field of : " << faItems_[i].meanFieldName() << "  with the type of " << SurfaceFieldType::typeName << nl;
                if (obrTgt_.found("mapped"+faItems_[i].meanFieldName()))
                {
                    obrTgt_.checkOut(*obrTgt_["mapped"+faItems_[i].meanFieldName()]);
                }
                addMappedFieldType<SurfaceFieldType,Type>(i);
                mapMeanSFieldType<SurfaceFieldType,Type>(i);
            }
        }
    }

}

template<class Type1, class Type2, class Type>
void Foam::functionObjects::mapAveField::mapPrime2MeanFieldType(const label fieldi)
{
    const Type1& mappedField = obrTgt_.lookupObject<Type1>("mapped"+faItems_[fieldi].meanFieldName());
    const Type2& fieldSqrMean = obr_.lookupObject<Type2>(faItems_[fieldi].fieldName()+"SqrMean");

    Type2& mappedPrimeField = const_cast<Type2&>
      (
            obrTgt_.lookupObject<Type2>("mapped"+faItems_[fieldi].prime2MeanFieldName())
      );

    cfdTorcfdInterpPtr->mapSrcToTgt(fieldSqrMean,plusEqOp<Type>(),mappedPrimeField);
    //mappedPrimeField = mappedPrimeField - sqr(mappedField);
    mappedPrimeField = mappedPrimeField - magSqr(mappedField);

}

template<class Type1, class Type2>
void Foam::functionObjects::mapAveField::mapPrime2MeanField()
{
    typedef GeometricField<Type1, fvPatchField, volMesh> VolFieldType1;
    typedef GeometricField<Type2, fvPatchField, volMesh> VolFieldType2;

    forAll(faItems_, i)
    {
        if (faItems_[i].prime2Mean())
        {
            const word& fieldName = faItems_[i].fieldName();
            if (obr_.foundObject<VolFieldType1>(fieldName))
            {
                mapPrime2MeanFieldType<VolFieldType1, VolFieldType2, Type2>(i);
            }
        }
    }
}

//========================= writing fields ===================================//

template<class Type>
void Foam::functionObjects::mapAveField::writeMappedFieldType
(
    const word& fieldName
) const
{

    if (obrTgt_.foundObject<Type>(fieldName))
    {
        Log << "    Write the mapped field :  " << fieldName << endl;
        const Type& f = obrTgt_.lookupObject<Type>(fieldName);
        f.write();
    }

}


template<class Type>
void Foam::functionObjects::mapAveField::writeMappedField() const
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> SurfaceFieldType;

    forAll(faItems_, fieldi)
    {
        if (faItems_[fieldi].mean())
        {
            const word& mappedFieldName = "mapped"+faItems_[fieldi].meanFieldName();
            writeMappedFieldType<VolFieldType>(mappedFieldName);
            writeMappedFieldType<SurfaceFieldType>(mappedFieldName);
        }
        if (faItems_[fieldi].prime2Mean())
        {
            const word& mappedFieldName = "mapped"+faItems_[fieldi].prime2MeanFieldName();
            writeMappedFieldType<VolFieldType>(mappedFieldName);
        }
    }
}

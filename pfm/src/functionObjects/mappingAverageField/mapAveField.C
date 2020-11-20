/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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
#include "mapAveField.H"
#include "fvCFD.H"
#include "fvMesh.H"
#include "meshToMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fieldAverageItem.H"
#include "autoPtr.H"
#include "meshSearch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(mapAveField, 0);
    addToRunTimeSelectionTable(functionObject, mapAveField, dictionary);
}
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::mapAveField::addMappedfields()
{

    Log << "    Initialising mapped fields" << endl;

    forAll(faItems_, i)
    {
        if (faItems_[i].mean())
        {
            if (obrTgt_.found("mapped"+faItems_[i].meanFieldName()))
            {
                obrTgt_.checkOut(*obrTgt_["mapped"+faItems_[i].meanFieldName()]);
            }
        }

        if (faItems_[i].prime2Mean())
        {
            if (obrTgt_.found("mapped"+faItems_[i].prime2MeanFieldName()))
            {
                obrTgt_.checkOut(*obrTgt_["mapped"+faItems_[i].prime2MeanFieldName()]);
            }
            if (obr_.found(faItems_[i].fieldName()+"SqrMean"))
            {
                obr_.checkOut(*obr_[faItems_[i].fieldName()+"SqrMean"]);
            }
        }
    }

    forAll(faItems_, fieldi)
    {

        addMappedField<scalar>(fieldi);
        addMappedField<vector>(fieldi);

    }

    forAll(faItems_, fieldi)
    {

        addMeanFieldSqr<scalar, scalar>(fieldi);
        addMeanFieldSqr<vector, scalar>(fieldi);
        addMappedPrimeField<scalar, scalar>(fieldi);
        addMappedPrimeField<vector, scalar>(fieldi);

    }

    mappingInitialised_=true;

    Log << endl;

}

void Foam::functionObjects::mapAveField::calculateAverageSqr()
{
    Log << "    calculating mean field-square" << endl;

    calcMeanFieldSqr<scalar, scalar>();
    calcMeanFieldSqr<vector, scalar>();

    Log << endl;
}

void Foam::functionObjects::mapAveField::mapAveragedField()
{
    Log << "    Mapping mean fields" << endl;

        mapMeanFields<scalar>();
        mapMeanFields<vector>();

        mapPrime2MeanField<scalar, scalar>();
        mapPrime2MeanField<vector, scalar>();

    Log << endl;

}

void Foam::functionObjects::mapAveField::writeMappedFields() const
{
    Log << "    Writing mapped fields" << endl;

        writeMappedField<scalar>();
        writeMappedField<vector>();

    Log << endl;
}

void Foam::functionObjects::mapAveField::writeMappedAverages()
{

    mapAveragedField();
    writeMappedFields();

}
// works only for scalars
void Foam::functionObjects::mapAveField::calcNewField()
{

    word newFieldName;
    wordList fieldsList_;
    scalarList factorsList_;

    IOdictionary mapFieldProperties
    (
        IOobject
        (
            "mapFieldsDict",
            rMesh_.time().system(),
            rMesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    mapFieldProperties.lookup("fieldName") >> newFieldName;
    mapFieldProperties.lookup("fieldsList") >> fieldsList_;
    mapFieldProperties.lookup("prefactorsList") >> factorsList_;

    const volScalarField& firstField = obrTgt_.lookupObject<volScalarField>(fieldsList_[0]);

    volScalarField newField
    (
        IOobject
        (
            newFieldName,
            rMesh_.time().timeName(),
            rMesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rMesh_,
        dimensioned<scalar>(newFieldName,firstField.dimensions(),0.0)
    );

    forAll(fieldsList_,i)
    {
        volScalarField addField = obrTgt_.lookupObject<volScalarField>(fieldsList_[i]);
        newField += factorsList_[i]*addField;
    }

    newField.write();

}

void Foam::functionObjects::mapAveField::findPais()
{


    const fvMesh& mesh_ = refCast<const fvMesh>(obr_);

    int nScells = mesh_.nCells();
    int nSfaces = mesh_.nFaces();

    int nTcells = rMesh_.nCells();
    int nTfaces = rMesh_.nFaces();

    Log << " n_source_cells = " << nScells << endl;
    Log << " n_source_faces = " << nSfaces << endl;

    Log << " n_target_cells = " << nTcells << endl;
    Log << " n_target_faces = " << nTfaces << endl;

    forAll(mesh_.cells(),cSrc)
    {

        Log << " source mesh cell ID = " << cSrc << endl;

        const cell& srcFaces = mesh_.cells()[cSrc];
        const vector& srcCentre = mesh_.cellCentres()[cSrc];

        label cTgt = rMesh_.findCell(srcCentre);
        Log << " target mesh cell ID = " << cTgt << endl;

        const cell& tgtFaces = rMesh_.cells()[cTgt];

        if(cTgt > -1)
        {
            Log << " source_faces = " << srcFaces << " target_faces = " << tgtFaces << endl;

            forAll(srcFaces,fi)
            {
                forAll(tgtFaces,fj)
                {
                    Log << " source face ID = " << srcFaces[fi] << " target face ID = " << tgtFaces[fj] << endl;
                    label sFace = srcFaces[fi];
                    label tFace = tgtFaces[fj];

                    const vector& ni= mesh_.faceAreas()[sFace] / mag(mesh_.faceAreas()[sFace]);
                    const vector& nj= rMesh_.faceAreas()[tFace] / mag(rMesh_.faceAreas()[tFace]);

                    if( mag( ni - nj) < SMALL || mag( ni + nj) < SMALL )
                    {

                        scalar sDist = ni & mesh_.faceCentres()[sFace];
                        scalar tDist = ni & rMesh_.faceCentres()[tFace];

                        label faceOwn = rMesh_.faceOwner()[tFace];

                        if( mag(sDist - tDist) < SMALL && (faceOwn == cTgt))
                        {
                            //add fi face to fj
                            Log << " source face : " << sFace << " matches with target face : " << tFace << endl;
                            labelPair facePair(sFace,tFace);
                            facePairs_.append(facePair);
                        }
                    }
                }
            }
        }

        else
        {
            FatalErrorInFunction
                << " No match for cell : " << cSrc
                << abort(FatalError);
        }
    }

    List<labelPairList> facePairList(Pstream::nProcs());
    facePairList[Pstream::myProcNo()] = facePairs_;
    Pstream::gatherList(facePairList);

    OFstream pairFile("facePairs");
    pairFile << facePairList;

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::mapAveField::mapAveField
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
  fieldAverage(name, runTime, dict),
  mappingInitialised_(false),
  addNewField_(false),
  Foam::fvMesh
    (
            IOobject
            (
                    "rMesh",
                    runTime.timeName(),
                    runTime,
                    IOobject::MUST_READ
            )
    ),

  obrTgt_
  (
      runTime.lookupObject<objectRegistry>
      (
          "rMesh"
      )
  ),
  rMesh_(refCast<const fvMesh>(obrTgt_)),

  facePairs_(0),

  cfdTorcfdInterpPtr(NULL)

{
    read(dict);
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::mapAveField::~mapAveField()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::mapAveField::read(const dictionary& dict)
{

    fieldAverage::read(dict);

    mappingInitialised_=false;

    const fvMesh& mesh_ = refCast<const fvMesh>(obr_);

    dict.lookup("addNewField") >> addNewField_;

    IOdictionary interpolationProperties
    (
        IOobject
        (
            "interpolationProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    const bool interpolatePatches(true);

    bool mapFlux;
    word mappingMethod;

    interpolationProperties.lookup("mapSurfaceFields") >> mapFlux;
    interpolationProperties.lookup("mappingMethod") >> mappingMethod;

    if(mappingMethod == "consistent")
    {
        if (rMesh_.bounds().overlaps(mesh_.bounds()))
        {

            cfdTorcfdInterpPtr =
            (
                new meshToMesh
                (

                    mesh_,  // source
                    rMesh_, // target
                    meshToMesh::interpolationMethodNames_.read
                    (
                        interpolationProperties.lookup("interpolationMethod")
                    ),
                    interpolatePatches // interpolate patches
                )
            );
         }

        else
        {
            FatalErrorInFunction
                << "regions " << rMesh_.name() << " and "
                << mesh_.name() <<  " do not intersect"
                << exit(FatalError);
        }
    }

    else if(mappingMethod == "cuttingMesh")
    {
        HashTable<word> patchMap;
        wordList cuttingPatches;

        IOdictionary mapFieldProperties
        (
            IOobject
            (
                "mapFieldsDict",
                rMesh_.time().system(),
                rMesh_,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
        );

        mapFieldProperties.lookup("patchMap") >> patchMap;
        mapFieldProperties.lookup("cuttingPatches") >>  cuttingPatches;


            cfdTorcfdInterpPtr =
            (
                new meshToMesh
                (

                    mesh_,  // source mesh
                    rMesh_, // target mesh
                    meshToMesh::interpolationMethodNames_.read
                    (
                        interpolationProperties.lookup("interpolationMethod")
                    ),
                    patchMap,
                    cuttingPatches
                )
            );
    }

    else
    {
        FatalErrorInFunction << "unknown mapping method " << mappingMethod << nl
            << " if the meshes are compatible choose consistent" << nl
            << " if the meshes are incompatible choose cuttingMesh" << nl
            << abort(FatalError);
    }

    if(mapFlux)
    {
        findPais();
    }


    return true;
}

bool Foam::functionObjects::mapAveField::execute()
{

    fieldAverage::execute();

    if (!mappingInitialised_)
    {
        addMappedfields();
    }

    calculateAverageSqr();

    return true;
}

bool Foam::functionObjects::mapAveField::write()
{


    writeMappedAverages();

    if(addNewField_)
    {
        calcNewField();
    }

    writeAveragingProperties();

    if (restartOnOutput_)
    {
        restart(); //must be comment out when fieldAverage::write() is written
        addMappedfields();
    }


    //fieldAverage::write();  //if sourceMesh mean fields are needed as well

    return true;
}


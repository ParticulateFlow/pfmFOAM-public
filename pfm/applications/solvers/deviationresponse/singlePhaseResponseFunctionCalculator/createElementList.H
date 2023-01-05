label locCellsOrFaces = 0;
const polyBoundaryMesh& boundaryMesh = mesh.boundaryMesh();
labelList faceIDperPatch(0,0);
labelList patchOwningFace(0,0);
if (internalSource)
{
    locCellsOrFaces = mesh.nCells();
}
else
{
    for (int i=0; i < patches.size(); i++)
    {
        label patchID = boundaryMesh.findPatchID(patches[i]);
        locCellsOrFaces += boundaryMesh[patchID].size();
    }
    faceIDperPatch.resize(locCellsOrFaces);
    patchOwningFace.resize(locCellsOrFaces);
    
    label fc = 0;
    for (int i=0; i < patches.size(); i++)
    {
        label patchID = boundaryMesh.findPatchID(patches[i]);

        forAll (boundaryMesh[patchID],faceI)
        {
            faceIDperPatch[fc] = faceI;
            patchOwningFace[fc] = patchID;
            fc++;
        }
    }
};

globalIndex globalNumbering(locCellsOrFaces);

labelList elementsLoc(locCellsOrFaces,0);
// fill elements list with GLOBAL IDs
for (int element = 0; element < locCellsOrFaces; element++)
{
    elementsLoc[element] = globalNumbering.toGlobal(element);
}


labelList numElementsLoc(Pstream::nProcs(),0);
numElementsLoc[Pstream::myProcNo()] = elementsLoc.size();
    
Pstream::gatherList(numElementsLoc);
Pstream::scatterList(numElementsLoc);

label totNumElements = 0;
for (int i = 0; i < numElementsLoc.size(); i++)
{
    totNumElements += numElementsLoc[i];
}

label accNumElements = 0;
for (int i = 0; i < Pstream::myProcNo(); i++)
{
    accNumElements += numElementsLoc[i];
}

labelList sourceElements(totNumElements,0);

for (int i = 0; i < elementsLoc.size(); i++)
{
    sourceElements[accNumElements + i] = elementsLoc[i];
}

Pstream::listCombineGather(sourceElements, plusEqOp<label>());
Pstream::listCombineScatter(sourceElements);

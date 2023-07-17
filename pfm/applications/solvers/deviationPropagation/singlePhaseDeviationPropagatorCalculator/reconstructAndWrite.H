// once all components are assembled, check for non-zero values
globalIndex globalCellNumbering(mesh.nCells());
labelList targetCells_loc(mesh.nCells(), -1); // *global* cell IDs on each processor
scalarList K_uu_loc(numCmpt * numCmpt * mesh.nCells(), 0.0);

scalar norm;
label numEntries = 0;
label iC = 0;
label cmpt = -1;
forAll(K_uu, cellI)
{
    norm = 0.0;
    forAll(components,c)
    {
        norm += magSqr(K_uu_allCmpt[c][cellI]);
    }
    if (Foam::sqrt(norm) * normalization > minKuu)
    {
        Receiver[cellI] = 1.0;
        targetCells_loc[numEntries] = globalCellNumbering.toGlobal(cellI);
        iC = 0;
        // OpenFOAM writes tensors as lists with components xx xy xz yx yy yz zx zy zz,
        // hence the inner loop goes over the elements of the list K_uu_allCmpt with fixed
        // component of the vectors K_uu computet in the main function
        forAll(components,c1)
        {
            cmpt = components[c1];
            forAll(components,c2)
            {
                K_uu_loc[numCmpt * numCmpt * numEntries + iC++] = K_uu_allCmpt[c2][cellI].component(cmpt);
            }
        }
        numEntries++;
    }
}
Receiver.write();

labelList numEntriesLoc(Pstream::nProcs(),0);
numEntriesLoc[Pstream::myProcNo()] = numEntries;

Pstream::gatherList(numEntriesLoc);
Pstream::scatterList(numEntriesLoc);

label totNumEntries = 0;
for (int i = 0; i < numEntriesLoc.size(); i++)
{
    totNumEntries += numEntriesLoc[i];
}

label accNumEntries = 0;
for (int i = 0; i < Pstream::myProcNo(); i++)
{
    accNumEntries += numEntriesLoc[i];
}

labelList targetCells_glob(totNumEntries, 0);
scalarList K_uu_glob(numCmpt * numCmpt * totNumEntries, 0.0);

for (int i = 0; i < numEntries; i++)
{
    targetCells_glob[accNumEntries + i] = targetCells_loc[i];
}

for (int i = 0; i < numCmpt*numCmpt*numEntries; i++)
{
    K_uu_glob[numCmpt*numCmpt*accNumEntries + i] = K_uu_loc[i];
}

Pstream::listCombineGather(targetCells_glob, plusEqOp<label>());
Pstream::listCombineGather(K_uu_glob, plusEqOp<scalar>());

if (Pstream::master())
{
    if (mode != integrated)
    {
        OS_C << sourceElements[sourceElement];
        for (int i = 0; i < totNumEntries; i++)
        {
            OS_C << " " << targetCells_glob[i];
        }
        OS_C << endl;

        OS_K << sourceElements[sourceElement];
        for (int i = 0; i < numCmpt*numCmpt*totNumEntries; i++)
        {
            OS_K << " " << K_uu_glob[i];
        }
        OS_K << endl;
    }
    else
    {
        for (int i = 0; i < totNumEntries; i++)
        {
            OS_K << i;
            for (int j = 0; j < numCmpt*numCmpt; j++)
            {
                OS_K << " " << K_uu_glob[i*numCmpt*numCmpt + j];
            }
            OS_K << endl;
        }
    }
}

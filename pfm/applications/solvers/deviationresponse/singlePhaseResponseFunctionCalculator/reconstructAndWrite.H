// once all components are assembled, check for non-zero values
globalIndex globalCellNumbering(mesh.nCells());
labelList targetCells_loc(mesh.nCells(), -1); // *global* cell IDs on each processor
scalarList X_uu_loc(numCmpt * numCmpt * mesh.nCells(), 0.0);

scalar norm;
label numEntries = 0;
label iC = 0;
label cmpt = -1;
forAll(X_uu, cellI)
{
    norm = 0.0;
    forAll(components,c)
    {
        norm += magSqr(X_uu_allCmpt[c][cellI]);
    }
    if (norm > minXuu2)
    {
        targetCells_loc[numEntries] = globalCellNumbering.toGlobal(cellI);
        iC = 0;
        forAll(components,c1)
        {
            forAll(components,c2)
            {
                cmpt = components[c2];
                X_uu_loc[numCmpt * numCmpt * numEntries + iC++] = X_uu_allCmpt[c1][cellI].component(cmpt);
            }
        }
        numEntries++;
    }
}

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
scalarList X_uu_glob(numCmpt * numCmpt * totNumEntries, 0.0);

for (int i = 0; i < numEntries; i++)
{
    targetCells_glob[accNumEntries + i] = targetCells_loc[i];
}

for (int i = 0; i < numCmpt*numCmpt*numEntries; i++)
{
    X_uu_glob[numCmpt*numCmpt*accNumEntries + i] = X_uu_loc[i];
}

Pstream::listCombineGather(targetCells_glob, plusEqOp<label>());
Pstream::listCombineGather(X_uu_glob, plusEqOp<scalar>());

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

        OS_X << sourceElements[sourceElement];
        for (int i = 0; i < numCmpt*numCmpt*totNumEntries; i++)
        {
            OS_X << " " << X_uu_glob[i];
        }
        OS_X << endl;
    }
    else
    {
        for (int i = 0; i < totNumEntries; i++)
        {
            OS_X << i;
            for (int j = 0; j < numCmpt*numCmpt; j++)
            {
                OS_X << " " << X_uu_glob[i*numCmpt*numCmpt + j];
            }
            OS_X << endl;
        }
    }
}

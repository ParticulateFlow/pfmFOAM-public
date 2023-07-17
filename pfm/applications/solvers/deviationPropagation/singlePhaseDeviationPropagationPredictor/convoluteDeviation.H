forAll(deltaUEvolved, cellI)
{
    labelList &senderCells = db.exportDeviationPropagators().senderCellIDs(nearestRefState,cellI);
    tensorList &K_uu = db.exportDeviationPropagators().Kuu_internal(nearestRefState,cellI);

    deltaUEvolved[cellI] = vector::zero;
    scalar norm = 0.0;
    forAll(K_uu, sender)
    {
        norm = K_uu[sender].xx()*K_uu[sender].xx() +
        K_uu[sender].xy()*K_uu[sender].xy() +
        K_uu[sender].xz()*K_uu[sender].xz() +
        K_uu[sender].yx()*K_uu[sender].yx() +
        K_uu[sender].yy()*K_uu[sender].yy() +
        K_uu[sender].yz()*K_uu[sender].yz() +
        K_uu[sender].zx()*K_uu[sender].zx() +
        K_uu[sender].zy()*K_uu[sender].zy() +
        K_uu[sender].zz()*K_uu[sender].zz();
        if (Foam::sqrt(norm) * mesh.V()[senderCells[sender]] > minKuu)
        {
            deltaUEvolved[cellI] += K_uu[sender] & deltaU[senderCells[sender]] * mesh.V()[senderCells[sender]];
        }
    }

    labelList &senderBoundaryFaces = db.exportDeviationPropagators().senderBoundaryFaceIDs(nearestRefState,cellI);
    tensorList &K_uu_boundary = db.exportDeviationPropagators().Kuu_boundary(nearestRefState,cellI);
    labelList &faceIDperPatch = db.faceIDperPatch();
    labelList &patchOwningFace = db.patchOwningFace();

    label patchID;
    label faceID;
    norm = 0.0;
    forAll(K_uu_boundary, sender)
    {
        norm = K_uu_boundary[sender].xx()*K_uu_boundary[sender].xx() +
        K_uu_boundary[sender].xy()*K_uu_boundary[sender].xy() +
        K_uu_boundary[sender].xz()*K_uu_boundary[sender].xz() +
        K_uu_boundary[sender].yx()*K_uu_boundary[sender].yx() +
        K_uu_boundary[sender].yy()*K_uu_boundary[sender].yy() +
        K_uu_boundary[sender].yz()*K_uu_boundary[sender].yz() +
        K_uu_boundary[sender].zx()*K_uu_boundary[sender].zx() +
        K_uu_boundary[sender].zy()*K_uu_boundary[sender].zy() +
        K_uu_boundary[sender].zz()*K_uu_boundary[sender].zz();

        patchID = patchOwningFace[senderBoundaryFaces[sender]];
        faceID = faceIDperPatch[senderBoundaryFaces[sender]];
        if (Foam::sqrt(norm) * mesh.magSf().boundaryField()[patchID][faceID] > minKuu)
        {
            deltaUEvolved[cellI] += K_uu_boundary[sender] & deltaU.boundaryField()[patchID][faceID] * mesh.magSf().boundaryField()[patchID][faceID];
        }
    }
}

forAll(deltaUEvolved, cellI)
{
    labelList &senderCells = db.responseF().senderCellIDs(nearestRefState,cellI);
    tensorList &X_uu = db.responseF().Xuu_internal(nearestRefState,cellI);

    deltaUEvolved[cellI] = vector::zero;
    forAll(X_uu, sender)
    {
        deltaUEvolved[cellI] += X_uu[sender] & deltaU[senderCells[sender]] * mesh.V()[senderCells[sender]];
    }

    labelList &senderBoundaryFaces = db.responseF().senderBoundaryFaceIDs(nearestRefState,cellI);
    tensorList &X_uu_boundary = db.responseF().Xuu_boundary(nearestRefState,cellI);
    labelList &faceIDperPatch = db.faceIDperPatch();
    labelList &patchOwningFace = db.patchOwningFace();

    label patchID;
    label faceID;
    forAll(X_uu_boundary, sender)
    {
        patchID = patchOwningFace[senderBoundaryFaces[sender]];
        faceID = faceIDperPatch[senderBoundaryFaces[sender]];
        deltaUEvolved[cellI] += X_uu_boundary[sender] & deltaU.boundaryField()[patchID][faceID] * mesh.magSf().boundaryField()[patchID][faceID];
    }
}

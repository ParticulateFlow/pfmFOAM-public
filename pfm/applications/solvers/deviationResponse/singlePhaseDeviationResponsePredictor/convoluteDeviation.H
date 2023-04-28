forAll(deltaUEvolved, cellI)
{
    labelList &senderCells = db.responseF().senderCellIDs(nearestRefState,cellI);
    tensorList &X_uu = db.responseF().Xuu_internal(nearestRefState,cellI);

    deltaUEvolved[cellI] = vector::zero;
    scalar norm = 0.0;
    forAll(X_uu, sender)
    {
        norm = X_uu[sender].xx()*X_uu[sender].xx() +
        X_uu[sender].xy()*X_uu[sender].xy() +
        X_uu[sender].xz()*X_uu[sender].xz() +
        X_uu[sender].yx()*X_uu[sender].yx() +
        X_uu[sender].yy()*X_uu[sender].yy() +
        X_uu[sender].yz()*X_uu[sender].yz() +
        X_uu[sender].zx()*X_uu[sender].zx() +
        X_uu[sender].zy()*X_uu[sender].zy() +
        X_uu[sender].zz()*X_uu[sender].zz();
        if (Foam::sqrt(norm) * mesh.V()[senderCells[sender]] > minXuu)
        {
            deltaUEvolved[cellI] += X_uu[sender] & deltaU[senderCells[sender]] * mesh.V()[senderCells[sender]];
        }
    }

    labelList &senderBoundaryFaces = db.responseF().senderBoundaryFaceIDs(nearestRefState,cellI);
    tensorList &X_uu_boundary = db.responseF().Xuu_boundary(nearestRefState,cellI);
    labelList &faceIDperPatch = db.faceIDperPatch();
    labelList &patchOwningFace = db.patchOwningFace();

    label patchID;
    label faceID;
    norm = 0.0;
    forAll(X_uu_boundary, sender)
    {
        norm = X_uu_boundary[sender].xx()*X_uu_boundary[sender].xx() +
        X_uu_boundary[sender].xy()*X_uu_boundary[sender].xy() +
        X_uu_boundary[sender].xz()*X_uu_boundary[sender].xz() +
        X_uu_boundary[sender].yx()*X_uu_boundary[sender].yx() +
        X_uu_boundary[sender].yy()*X_uu_boundary[sender].yy() +
        X_uu_boundary[sender].yz()*X_uu_boundary[sender].yz() +
        X_uu_boundary[sender].zx()*X_uu_boundary[sender].zx() +
        X_uu_boundary[sender].zy()*X_uu_boundary[sender].zy() +
        X_uu_boundary[sender].zz()*X_uu_boundary[sender].zz();

        patchID = patchOwningFace[senderBoundaryFaces[sender]];
        faceID = faceIDperPatch[senderBoundaryFaces[sender]];
        if (Foam::sqrt(norm) * mesh.magSf().boundaryField()[patchID][faceID] > minXuu)
        {
            deltaUEvolved[cellI] += X_uu_boundary[sender] & deltaU.boundaryField()[patchID][faceID] * mesh.magSf().boundaryField()[patchID][faceID];
        }
    }
}

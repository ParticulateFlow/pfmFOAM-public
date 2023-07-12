#!/bin/sh

sed -i -e 's/volScalarField/labelList/g' 0/cellDist
sed -i -e 's/"0"/"constant"/g' 0/cellDist

sed -i -e 's/volScalarField/labelList/g' 0/rMesh/cellDist
sed -i -e 's/"0/"constant/g' 0/rMesh/cellDist

sed -i -e '/dimensions/d' 0/cellDist
sed -i -e '/dimensions/d' 0/rMesh/cellDist

sed -i -e '/internalField/,+1d' 0/cellDist
sed -i -e '/internalField/,+1d' 0/rMesh/cellDist

sed -i -e '/boundaryField/,$d' 0/cellDist
sed -i -e '/boundaryField/,$d' 0/rMesh/cellDist

cp 0/cellDist constant/
cp 0/rMesh/cellDist constant/rMesh/

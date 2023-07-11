#!/bin/bash
cd referenceDataGeneration
cp -r orig.0 0
blockMesh > blockMesh.log 2>&1
pisoFoam > pisoFoam.log 2>&1
mkdir data_referenceTimeSeries
mv [1-9]* data_referenceTimeSeries
cd data_referenceTimeSeries
rm -r */uniform
rm */phi
cd ../..

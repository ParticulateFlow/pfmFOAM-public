#!/bin/bash
cd validationData
cd uin_120
cp -r orig.0 0
cp system/controlDict_init system/controlDict
blockMesh > blockMesh.log 2>&1
pisoFoam > pisoFoam_init.log 2>&1

mkdir data_uin_120
cp system/controlDict_valid_multistep system/controlDict
pisoFoam > pisoFoam_valid_multistep.log 2>&1

mv [1-2]* data_uin_120/
cp -r data_uin_120/181 .

cp system/controlDict_valid_long system/controlDict
pisoFoam > pisoFoam_valid_long.log 2>&1
mv 2* data_uin_120/
mv postProcessing data_uin_120/
rm -r data_uin_120/[1-9]*/uniform
rm data_uin_120/[1-9]*/phi

cd ../..

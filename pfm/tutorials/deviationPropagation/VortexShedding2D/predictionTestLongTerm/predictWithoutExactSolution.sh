#!/bin/bash
path="../validationData/uin_120/data_uin_120"
t1=181
t2=281

mkdir log

rm -r 0
mkdir 0
echo "copying from $path/$t1 to 0"
cp $path/$t1/U 0/
# get p file to correct solution if necessary
cp $path/$t1/p 0/
rm distances

# run simulation
singlePhaseDeviationPropagationPredictor > log/run_$t1.log #2>&1

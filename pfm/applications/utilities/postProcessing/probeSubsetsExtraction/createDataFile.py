#!/usr/bin/env python3
import os
import sys

#################################################
# user-defined input
nLabels = 2

probeFieldsScalar = ["p"]
probeFieldsVector = ["U"]
#################################################

numProbeFieldsScalar = len(probeFieldsScalar)
numProbeFieldsVector = len(probeFieldsVector)

labels=[]

# get sample times
with open('sampleTimes', 'r') as f:
    sampleTimes = f.readlines()
numSampleTimes = len(sampleTimes)

# get extracted probe times
with open(os.path.join(sampleTimes[0].strip(),"pointEvol_"+probeFieldsScalar[0]+"_cellSet_0"), 'r') as f:
    counter = 0
    for line in f:
        if line.startswith('#'): 
            continue
        counter += 1
        if (counter>nLabels): break
        labels.append(line[0])


# get cell sets
with open('cellSetsExtensions') as f:
    numCellSets = 0
    for line in f:
        if line.startswith('#'): 
            continue
        numCellSets = numCellSets + 1;

db = open("dataBase.csv","w+")

# first write header line
with open('distData_header') as f:
    db.write(f.read());
for fieldI in probeFieldsScalar:
    for i in labels:
        db.write(", " + fieldI + "_evolved_" + i)
for fieldI in probeFieldsVector:
    for i in labels:
        db.write(", " + fieldI + "X_evolved_" + i)
for fieldI in probeFieldsVector:
    for i in labels:
        db.write(", " + fieldI + "Y_evolved_" + i)
for fieldI in probeFieldsVector:
    for i in labels:
        db.write(", " + fieldI + "Z_evolved_" + i)
db.write("\n");

# now insert data
for sampleTime in sampleTimes:
    for cellSetIndex in range(numCellSets):
        with open(os.path.join(sampleTime.strip(),"distData_cellSet_"+str(cellSetIndex)+".csv"), 'r') as fDist:
            db.write(fDist.read())
        for fieldI in probeFieldsScalar:
            with open(os.path.join(sampleTime.strip(),"pointEvol_"+fieldI+"_cellSet_"+str(cellSetIndex)), 'r') as fPointScalar:
                counter = 0
                for line in fPointScalar:
                    if line.startswith('#'): 
                        continue
                    counter += 1
                    if (counter>nLabels): break
                    db.write(", "+line[1].rstrip("\n"))
            if (counter<nLabels): sys.exit('Could not find enough probe values for sampleTime '+sampleTime)

        for fieldI in probeFieldsVector:
            with open(os.path.join(sampleTime.strip(),"pointEvol_"+fieldI+"X_cellSet_"+str(cellSetIndex)), 'r') as fPointVector:
                counter = 0
                for line in fPointVector:
                    if line.startswith('#'): 
                        continue
                    counter += 1
                    if (counter>nLabels): break
                    db.write(", "+line[1].rstrip("\n"))
            with open(os.path.join(sampleTime.strip(),"pointEvol_"+fieldI+"Y_cellSet_"+str(cellSetIndex)), 'r') as fPointVector:
                counter = 0
                for line in fPointVector:
                    if line.startswith('#'): 
                        continue
                    counter += 1
                    if (counter>nLabels): break
                    db.write(", "+line[1].rstrip("\n"))
            with open(os.path.join(sampleTime.strip(),"pointEvol_"+fieldI+"Z_cellSet_"+str(cellSetIndex)), 'r') as fPointVector:
                counter = 0
                for line in fPointVector:
                    if line.startswith('#'): 
                        continue
                    counter += 1
                    if (counter>nLabels): break
                    db.write(", "+line[1].rstrip("\n"))

        db.write("\n")

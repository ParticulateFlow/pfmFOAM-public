#!/usr/bin/env python3
import os
import sys


nTotSteps = 3
stepSize = 2

acc = 0.05

probesPath = "database/postProcessing/probes1/0/"

probeFieldsScalar = ["p"]
numProbeFieldsScalar = len(probeFieldsScalar)

probeFieldsVector = ["U"]
numProbeFieldsVector = len(probeFieldsVector)

# probe numbers starting from 0
probePointIndices = [1]

def searchProbeList(reftime,times,values):
    probevaluesubset = []
    numTimes = len(times)
    for timeIndex in range(numTimes):
        if (abs(times[timeIndex]-float(reftime)) < 1e-8):
            for i in range(timeIndex,timeIndex+nTotSteps,stepSize):
                if (i>=len(values)): sys.exit('Could not find enough probe values for reftime '+reftime)
                probevaluesubset.append(values[i])
            return probevaluesubset
    return []
        
def getCellSetIndex(f,probeIndex):
    counter = 0
    for line in f:
        if counter == probeIndex:
            row = line.split()
            print(line)
            break
        counter = counter + 1
    x = float(row[3].lstrip("("))
    y = float(row[4])
    z = float(row[5].rstrip(")"))
    counter = 0
    with open('cellSetsExtensions') as f2:
        for line in f2:
            if line.startswith('#'): 
                continue
            row = line.split()
            xmin = float(row[1])
            xmax = float(row[2])
            ymin = float(row[3])
            ymax = float(row[4])
            zmin = float(row[5])
            zmax = float(row[6])

            dx = (x - 0.5 * (xmin + xmax)) / (xmax - xmin)
            dy = (y - 0.5 * (ymin + ymax)) / (ymax - ymin)
            dz = (z - 0.5 * (zmin + zmax)) / (zmax - zmin)

            if (abs(dx) < acc and abs(dy) < acc and abs(dz) < acc):
                return counter

            counter = counter + 1

    sys.exit('Could not find cellSet with center close to probing point.')
    


with open('sampleTimes', 'r') as f:
    sampleTimes = f.readlines()
numSampleTimes = len(sampleTimes)

# scalar fields
for probePointIndex in probePointIndices:
    for fieldI in probeFieldsScalar:
        t = []
        probevalue = []
        filename=probesPath+fieldI
        if not os.path.exists(filename): break
        with open(filename) as f:
            cellSetIndex = getCellSetIndex(f,probePointIndex)
            for line in f:
	        if line.startswith('#'): 
                    continue
                row = line.split()
                t.append(float(row[0]))
    	        probevalue.append(float(row[probePointIndex+1]))

        for sampleTime in sampleTimes:
            probesubset = []
            probesubset = searchProbeList(sampleTime,t,probevalue)
            with open(os.path.join(sampleTime.strip(),"pointEvol_"+fieldI+"_cellSet_"+str(cellSetIndex)), 'w') as f:
                for item in range(len(probesubset)):
                    f.write("%d %s\n" % (item*stepSize,probesubset[item]))

# vector fields
for probePointIndex in probePointIndices:
    for fieldI in probeFieldsVector:
        t = []
        probevalueX = []
        probevalueY = []
        probevalueZ = []
        filename=probesPath+fieldI
        if not os.path.exists(filename): break
        with open(filename) as f:
            cellSetIndex = getCellSetIndex(f,probePointIndex)
            for line in f:
	        if line.startswith('#'): 
                    continue
                row = line.split()
                t.append(float(row[0]))
    	        probevalueX.append(float(row[3*probePointIndex+1].lstrip("(")))
    	        probevalueY.append(float(row[3*probePointIndex+2]))
    	        probevalueZ.append(float(row[3*probePointIndex+3].rstrip(")")))

        for sampleTime in sampleTimes:
            probesubset = []
            probesubset = searchProbeList(sampleTime,t,probevalueX)
            with open(os.path.join(sampleTime.strip(),"pointEvol_"+fieldI+"X_cellSet_"+str(cellSetIndex)), 'w') as f:
                for item in probesubset:
                    f.write("%s\n" % item)

            probesubset = []
            probesubset = searchProbeList(sampleTime,t,probevalueY)
            with open(os.path.join(sampleTime.strip(),"pointEvol_"+fieldI+"Y_cellSet_"+str(cellSetIndex)), 'w') as f:
                for item in probesubset:
                    f.write("%s\n" % item)

            probesubset = []
            probesubset = searchProbeList(sampleTime,t,probevalueZ)
            with open(os.path.join(sampleTime.strip(),"pointEvol_"+fieldI+"Z_cellSet_"+str(cellSetIndex)), 'w') as f:
                for item in probesubset:
                    f.write("%s\n" % item)

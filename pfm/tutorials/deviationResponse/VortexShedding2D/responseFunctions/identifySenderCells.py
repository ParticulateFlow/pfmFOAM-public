#!/usr/bin/python
import os
import numpy as np
from numpy import random
from bisect import bisect_left ,bisect
import math
import re


import os.path
from os import path
from importlib.machinery import SourceFileLoader
if path.exists("/home/k3b02/k3b1677/"):
  p = "/home/k3b02/k3b1677/"
else:
  p = "/home/lichtenegger/"
readDataFile = SourceFileLoader("readDataFile",p+"scripts/dataAnalysis/readDataFile.py").load_module()
dataAnalysis = SourceFileLoader("dataAnalysis",p+"scripts/dataAnalysis/dataAnalysis.py").load_module()

mode = 'internal' # 'boundary'
filename = 'testtargetCells_'+mode
field = 'testX_uu_'
filenameX = field+mode

components = [0,1]
num_cmpts = len(components)
num_cmpts2 = num_cmpts * num_cmpts

senders = []
targets = []
X = []


with open(filename) as f:
    for line in f:
        if line.startswith('#'):
            continue
        row = re.split(' ',line.lstrip())
        senders.append(int(row[0]))
        t = list(map(int,row[1:]))
        targets.append(t)

with open(filenameX) as f:
    for line in f:
        if line.startswith('#'):
            continue
        row = re.split(' ',line.lstrip())
#        senders.append(int(row[0]))
        x = list(map(int,row[1:]))
        X.append(x)


# find minimum and maximum target indices
i0 = 1000000000
i1 = 0
for i in range(len(targets)):
    if targets[i][-1] > i1:
        i1 = targets[i][-1]
    if targets[i][0] < i0:
        i0 = targets[i][0]

senders_by_targets = np.empty((i1-i0+1, 0)).tolist()
X_by_targets = np.empty((i1-i0+1, 0)).tolist()

for i in range(i0,i1+1):
    for j in range(len(targets)):
        left_index = bisect_left(targets[j], i)
        if left_index!=bisect(targets[j], i):
            senders_by_targets[i-i0].append(senders[j])
            for k in range(num_cmpts2):
                X_by_targets[i-i0].append(X[j][left_index * num_cmpts2 + k])

f = open('sendersCells_'+mode,'w')
f.write("%d\n(\n" % int(i1-i0+1))
for i in range(i0,i1+1):
    f.write("%d(" % len(senders_by_targets[i-i0]))
    for j in range(len(senders_by_targets[i-i0])):
        f.write("%d" % senders_by_targets[i-i0][j])
        if j < len(senders_by_targets[i-i0])-1:
            f.write(" ")
    f.write(")\n")
f.write(")")
f.close()

f = open('receiverCells_'+mode,'w')
f.write("%d\n(" % int(i1-i0+1))
for i in range(i0,i1+1):
    f.write("%d" % i)
    if i < i1:
        f.write(" ")
f.write(")")
f.close()

f = open(field+'by_targets_'+mode,'w')
f.write("%d\n(\n" % int(i1-i0+1))
for i in range(i0,i1+1):
    f.write("%d(" % len(senders_by_targets[i-i0]))
    for j in range(len(senders_by_targets[i-i0])):
        Xlist = np.zeros(9)
        counter = 0
        for c1 in range(num_cmpts):
            for c2 in range(num_cmpts):
                cc1 = components[c1]
                cc2 = components[c2]
                Xlist[cc1 * 3 + cc2] = X_by_targets[i-i0][j * num_cmpts2 + counter]
                counter += 1
        f.write("(")
        for k in range(9):
            f.write("%f" % Xlist[k])
            if k < 8:
                f.write(" ")
        f.write(")")
        if j < len(senders_by_targets[i-i0]) - 1:
            f.write(" ")
    f.write(")\n")
f.write(")")
f.close()

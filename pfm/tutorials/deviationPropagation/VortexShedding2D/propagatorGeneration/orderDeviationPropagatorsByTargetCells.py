#!/usr/bin/python
import os
import numpy as np
from numpy import random
from bisect import bisect_left ,bisect
import math
import re
import sys
import os.path
from os import path

components = [0,1]
dataBaseName = 'dataBase'

def orderResponseFunctionsByTargetCells(path,mode,components):
    num_cmpts = len(components)
    num_cmpts2 = num_cmpts * num_cmpts
    filenameC = path+'targetCells_'+mode
    field = 'K_uu_'
    filenameX = path+field+'bySenders_'+mode
    senders = []
    targets = []
    X = []

    with open(filenameC) as f:
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
            x = list(map(float,row[1:]))
            X.append(x)

# find minimum and maximum target indices
# OPTION A: use full range of global cell indices by counting lines in the response to internal deviations
    i0 = 0
    i1 = 0
    with open(path+field+'bySenders_internal') as f:
        for line in f:
            if line.startswith('#'):
                continue
            i1 += 1
    i1 -= 1

# OPTION B: look for minimum and maximum targets
#           gives shorter list but does not necessarily contain all cells

#    i0 = 1000000000
#    i1 = 0
#    for i in range(len(targets)):
#        if len(targets[i]) == 0:
#            continue
#        if targets[i][-1] > i1:
#            i1 = targets[i][-1]
#        if targets[i][0] < i0:
#            i0 = targets[i][0]

    senders_by_targets = np.empty((i1-i0+1, 0)).tolist()
    X_by_targets = np.empty((i1-i0+1, 0)).tolist()

    for i in range(i0,i1+1):
        for j in range(len(targets)):
            left_index = bisect_left(targets[j], i)
            if left_index!=bisect(targets[j], i):
                senders_by_targets[i-i0].append(senders[j])
                for k in range(num_cmpts2):
                    X_by_targets[i-i0].append(X[j][left_index * num_cmpts2 + k])

    filename = path+'senderCells'
    if mode == 'boundary':
        filename = path+'senderFaces'
    f = open(filename,'w')
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

    f = open(path+'receiverCells_'+mode,'w')
    f.write("%d\n(" % int(i1-i0+1))
    for i in range(i0,i1+1):
        f.write("%d" % i)
        if i < i1:
            f.write(" ")
    f.write(")")
    f.close()

    f = open(path+field+mode+'.txt','w')
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

for refState in os.listdir(dataBaseName):
    path=dataBaseName+'/'+refState+'/'
    print(path)
    orderResponseFunctionsByTargetCells(path,'internal',components)
    orderResponseFunctionsByTargetCells(path,'boundary',components)

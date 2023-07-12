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

def formatIntegratedResponseFunctions(path,components):
    num_cmpts = len(components)
    num_cmpts2 = num_cmpts * num_cmpts
    field = 'K_uu_'
    mode = 'integrated'
    filenameX = path+field+'bySenders_'+mode
    senders = []
    targets = []
    X = []

    with open(filenameX) as f:
        for line in f:
            if line.startswith('#'):
                continue
            row = re.split(' ',line.lstrip())
            x = list(map(float,row[1:]))
            X.append(x)

    f = open(path+field+mode+'.txt','w')
    f.write("%d\n(" % len(X))

    for i in range(len(X)):
        Xlist = np.zeros(9)
        counter = 0
        f.write('(')
        for c1 in range(num_cmpts):
            for c2 in range(num_cmpts):
                cc1 = components[c1]
                cc2 = components[c2]
                Xlist[cc1 * 3 + cc2] = X[i][counter]
                counter += 1
        for j in range(9):
            f.write('%f' % Xlist[j])
            if j < 8:
                f.write(' ')
        f.write(')')

    f.write(')')
    f.close()

for refState in os.listdir(dataBaseName):
    path=dataBaseName+'/'+refState+'/'
    print(path)
    formatIntegratedResponseFunctions(path,components)

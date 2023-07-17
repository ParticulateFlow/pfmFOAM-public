#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import re
from matplotlib import rc
from matplotlib.ticker import AutoMinorLocator

import os.path
from os import path

def readDataFile(filename,col,split_at=r'\s+',ignore='#',ignoreleadinglines=0):
    data_ = []
    counter = 0
    with open(filename) as f:
        for line in f:
            if counter < ignoreleadinglines:
                counter += 1
                continue
            if line.startswith(ignore):
                continue
            row = re.split(split_at,line.lstrip().rstrip())
            data_.append(float(row[col]))

    data = np.array(data_)
    return data

rc('font',**{'family':'serif','serif':['Times'],'size':5})
rc('text', usetex=True)

# number of steps in each sample
nsteps = 20


err_pred_list_uin120 = readDataFile('predictionTestMultiStep/distancesList',3)
err_pred_uin120 = np.zeros(nsteps)
std_pred_uin120 = np.zeros(nsteps)

for i in range(nsteps):
   subarray = err_pred_list_uin120[i::nsteps]
   err_pred_uin120[i] = np.mean(subarray)
   std_pred_uin120[i] = np.std(subarray)


times = np.arange(1,21,1)

fig = plt.figure()

plt.xlim([0,20])
plt.ylim([0,0.05])

ax1 = fig.add_subplot(111)

ax1.set_xlabel(r'$\displaystyle n_{\textrm{steps}}$',size=9)
ax1.set_ylabel(r'$\displaystyle E_{\textrm{pred}}$',size=9)
ax1.xaxis.set_minor_locator(AutoMinorLocator(2))
ax1.yaxis.set_minor_locator(AutoMinorLocator(2))

plt.xticks(np.arange(0, 20.001, 5.0))
plt.yticks(np.arange(0, 0.05001, 0.025))

ax1.errorbar(times,err_pred_uin120, yerr=std_pred_uin120, c='r',label=r'$\displaystyle u_{\textrm{in}} = 1.2$',fmt='s',markersize=1.0,elinewidth=0.5)

leg = ax1.legend(loc='upper left',frameon=False,ncol=3,handletextpad=0.5,)

plt.subplots_adjust(top=0.95,bottom=0.2,left=0.2,right=0.84)

fig.set_size_inches(3, 2.1)
plt.savefig("figures/multiStepPrediction.pdf", dpi=1000)
plt.savefig("figures/multiStepPrediction.png", dpi=1000)

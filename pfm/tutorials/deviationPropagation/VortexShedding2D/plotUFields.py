#!/usr/bin/python3

import numpy as np
import vtk
import matplotlib.pyplot as P
from matplotlib.ticker import StrMethodFormatter
from matplotlib import rc
from matplotlib import rcParams
import os
from os import path
import re

def loadVTK(filename,dim1,dim2,dimData):
    if not os.path.exists(filename): return None
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(filename)
    reader.ReadAllScalarsOn()
    reader.ReadAllVectorsOn()
    reader.Update()

    data = reader.GetOutput()
    cells = data.GetPolys()
    triangles = cells.GetData()
    points = data.GetPoints()
    point_data = data.GetPointData()
    Edata = point_data.GetArray(0)

    ntri = int(triangles.GetNumberOfTuples()/4)
    npts = points.GetNumberOfPoints()
    nvls = Edata.GetNumberOfTuples()

    tri = np.zeros((ntri, 3))
    x = np.zeros(npts)
    y = np.zeros(npts)
    data = np.zeros(nvls)

    for i in range(0, ntri):
        tri[i, 0] = triangles.GetTuple(4*i + 1)[0]
        tri[i, 1] = triangles.GetTuple(4*i + 2)[0]
        tri[i, 2] = triangles.GetTuple(4*i + 3)[0]

    for i in range(points.GetNumberOfPoints()):
        pt = points.GetPoint(i)
        x[i] = pt[dim1]
        y[i] = pt[dim2]

    for i in range(0, Edata.GetNumberOfTuples()):
        E = Edata.GetTuple(i)
        data[i] = E[dimData]
 

    return (x, y, tri, data)


rc('font',**{'family':'serif','serif':['Times'],'size':7})
rc('text', usetex=True)
rcParams['text.latex.preamble'] = [r'\usepackage{bm}']

if not os.path.exists("figures"):
     os.makedirs("figures")

P.clf()
P.figure(figsize=(2.65,1.86),dpi=1000)

P.xlabel('$x$ [m]')
xtics=np.linspace(-20, 30, 6)
P.gca().set_xticks(xtics)

P.ylabel('$y$ [m]')
ytics=np.linspace(-10, 10.0, 3)
P.gca().set_yticks(ytics)

P.xlim(-10, 20.001)
P.ylim(-10, 10.001)

P.subplots_adjust(left=0.16, right=0.98, top=0.90, bottom=0.23)

#################
# validation data
#################
x, y, tri, ux = loadVTK('validationData/uin_120/data_uin_120/postProcessing/cuttingPlane/281/UMean_zNormal.vtk',0,1,0)
x, y, tri, uy = loadVTK('validationData/uin_120/data_uin_120/postProcessing/cuttingPlane/281/UMean_zNormal.vtk',0,1,1)
uabs = np.sqrt(np.square(ux)+np.square(uy))

llevels = np.linspace(0,1.5, 10)
clevels = np.linspace(0, 1.5, 13)
cnt = P.tricontourf(x, y, tri, uabs, clevels, cmap=P.cm.coolwarm, extend="max",extendfrac=0.01)

v = np.linspace(0.0, 1.5, 2, endpoint=True)
cbar = P.colorbar(ticks=v,extendrect=True,format='%.1f')
cbar.set_label(r'$\displaystyle |\bm{u}| [\textrm{m/s}]$',labelpad=-5)

P.savefig('figures/UMeanExact.png', dpi=1000)
P.savefig('figures/UMeanExact.pdf', dpi=1000)

##################################
# deviation-propagation prediction
##################################
x, y, tri, ux = loadVTK('predictionTestLongTerm/postProcessing/cuttingPlane/100/UMean_zNormal.vtk',0,1,0)
x, y, tri, uy = loadVTK('predictionTestLongTerm/postProcessing/cuttingPlane/100/UMean_zNormal.vtk',0,1,1)
uabs = np.sqrt(np.square(ux)+np.square(uy))

cnt = P.tricontourf(x, y, tri, uabs, clevels, cmap=P.cm.coolwarm, extend="max",extendfrac=0.01)

P.savefig('figures/UMeanDeviationPropagation.png', dpi=1000)
P.savefig('figures/UMeanDeviationPropagation.pdf', dpi=1000)

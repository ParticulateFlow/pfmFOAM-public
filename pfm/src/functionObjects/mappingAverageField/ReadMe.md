# mappingAveragedField

## Utility Intoduction
mapAveField is developed from the fieldAverage functionObject of OpenFOAM-4.x in order to map the averaged fields during the simulation from the fine grid to a coarse mesh.

For the surface fields this functionObeject finds the face pairs of the fine and coarse grids which are exactly on the same plane and belong to the corresponding cells. In mesh generation it is necessary to make sure that fine and coarse grids have the same nodes and faces.

Additionaly, for running in parallel it is required to make sure all of the corresponding cells in fine and coarse grids are decomposed into the same processor. Therefore decomposing should be performed manually. In order to do the manual decomposition we specify the cell sets with setFields.

## Compile the functionObject
For compiling the utility OpenFOAM-6 should be installed.

## Tutorial Case
To see how it works there is a cylinder_mapAveField tutorial case. Also, manual decomposing is included.

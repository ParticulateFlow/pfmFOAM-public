# twoPhaseEulerTurbFoam
OpenFOAM implementation of ADM and the SATFM model for the coarse grid simulation of gas-solid flows.
"twoPhaseEulerTurbFoam" is based on the OpenFOAM solver "twoPhaseEulerFoam".

## OpenFOAM compatibilty
twoPhaseEulerTurbFoam is currently compatible with the OpenFOAM-6 version including the commit
'thermophysicalModels: Added laminar thermal diffusivity for energy, alphahe
Needed for laminar transport of he (h or e)
Latest commit 1a0c91b on 5 Aug 2018'

## Compiling Instructions
Before compiling "twoPhaseEulerTurbFoam" the user lib "customLESfilters" has to be compiled.
Thus, execute 
```
wmake libso
```
in
```
pfm/src/TurbulenceModels/turbulenceModels
```

Finally, use 
```
./Allclean
./Allwmake
```
to compile the solver.

## Tutorials
For both cases (ADM and SATFM) a simple fluidized bed tutorial is provided (compare with Schneiderbauer & Saeedipour, 2019).

## References
* Schneiderbauer, S. (2017). A spatially-averaged two-fluid model for dense large-scale gas-solid flows. AIChE Journal, 63(8), 3544–3562.
* Schneiderbauer, S. (2018). Validation study on spatially averaged two-fluid model for gas-solid flows. II: Application to risers and bubbling fluidized beds. AIChE Journal, 64(5), 1606–1617. 
* Schneiderbauer, S. (2018). Validation study on spatially averaged two-fluid model for gas-solid flows: I. A-priori analysis of wall bounded flows. AIChE Journal, 64(5), 1591–1605. 
* Schneiderbauer, S., & Saeedipour, M. (2018). Approximate deconvolution model for the simulation of turbulent gas-solid flows: An a-priori analysis. Physics of Fluids, 30(2), 023301.
* Schneiderbauer, S., & Saeedipour, M. (2019). Numerical simulation of turbulent gas-solid flow using an approximate deconvolution model. International Journal of Multiphase Flow, 114, 287–302.

# twoPhaseEulerTurbFoam
OpenFOAM implementation of ADM and the SATFM model for the coarse grid simulation of gas-solid flows.
"twoPhaseEulerTurbFoam" is based on the OpenFOAM solver "twoPhaseEulerFoam".

## OpenFOAM compatibilty
twoPhaseEulerTurbFoam is currently compatible with the OpenFOAM-6 version including the commit
'thermophysicalModels: Added laminar thermal diffusivity for energy, alphahe
Needed for laminar transport of he (h or e)
Latest commit 1a0c91b on 5 Aug 2018'

## SATFM model for heat transfer
In addition to the SATFM correction to the momentum balance, a correction to the thermal energy balance in form of a drift temperature is provided.

Please note that due to the current implementation, phase 1 should be the dispersed particle-phase and phase 2 the gas-phase if the energy equation is solved in the SATFM model.

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
For both cases (ADM and SATFM) a simple fluidized bed tutorial is provided (compare with Schneiderbauer & Saeedipour, 2019). In addition, tutorials are added, which showcase the predictions of the thermal energy balance eautions for both, kinetic-theory Two-Fluid Model and SATFM (compare with Rauchenzauner & Schneiderbauer, Int. J. Heat Mass Transf., 2022).

## References
* Schneiderbauer, S. (2017). A spatially-averaged two-fluid model for dense large-scale gas-solid flows. AIChE Journal, 63(8), 3544–3562.
* Schneiderbauer, S. (2018). Validation study on spatially averaged two-fluid model for gas-solid flows. II: Application to risers and bubbling fluidized beds. AIChE Journal, 64(5), 1606–1617. 
* Schneiderbauer, S. (2018). Validation study on spatially averaged two-fluid model for gas-solid flows: I. A-priori analysis of wall bounded flows. AIChE Journal, 64(5), 1591–1605. 
* Schneiderbauer, S., & Saeedipour, M. (2018). Approximate deconvolution model for the simulation of turbulent gas-solid flows: An a-priori analysis. Physics of Fluids, 30(2), 023301.
* Schneiderbauer, S., & Saeedipour, M. (2019). Numerical simulation of turbulent gas-solid flow using an approximate deconvolution model. International Journal of Multiphase Flow, 114, 287–302.
* Rauchenzauner, S. and Schneiderbauer, S. (2020). A Dynamic Anisotropic Spatially-Averaged Two-Fluid Model for Moderately Dense Gas-Particle Flows. International Journal of Multiphase Flow, 126, 103237.
* Rauchenzauner, S. and Schneiderbauer, S. (2022). A dynamic multiphase turbulence model for coarse-grid simulations of fluidized gas-particle suspensions. Chemical Engineering Science, 247, 117104.
* Rauchenzauner, S. and Schneiderbauer, S. (2020). A Dynamic Spatially-Averaged Two-Fluid Model for Heat Transport in Moderately Dense Gas-Particle Flows. Physics of Fluids, 32, 063307.
* Rauchenzauner, S. and Schneiderbauer, S. (2022). Validation study of a Spatially-Averaged Two-Fluid Model for heat transport in gas-particle flows. International Journal of Heat and Mass Transfer, 198, 123382.

### The case setup for data generation is adapted from http://www.wolfdynamics.com/tutorials.html

Requirements: successfully compiled pfmFOAM; python3 installation

1. Generation of reference data: execute generateReferenceData.sh;
   the resulting time series is moved to referenceData/data_referenceTimeSeries
2. Generation of validation data (uniform inlet velocity 1.2 m/s):
   execute generateValidationData.sh;
   validation data is moved to validationData/data_uin_120
3. Generation of propagators: execute generatePropagator.sh;
   for each time in propagatorGeneration/listOfReferenceStates, the corresponding propagator is computed;
   this can take several days - to speed up the process, several copies
   propagatorGeneration0, propagatorGeneration1 etc. of propagatorGeneration
   may be created and the elements of listOfReferenceStates distributed among them;
   computations can be carried out in parallel and the resulting propagators found in
   propagatorGeneration0/dataBase, propagatorGeneration1/dataBase copied into one folder propagatorGeneration/dataBase
4. If propagators have been successfully computed and converted into binary format
   (in addition to files X_uu_internal.txt, there are also files X_uu_internal etc.),
   intermediate files and propagators in ASCII format may be deleted to free disk space.
   To do so, execute cleanDatabase.sh
5. Predictions for up to 20 steps for a case with increased inlet velocity u_in = 1.2 m/s:
   execute predictMultiSteps.sh;
   resulting data can be plotted with "python3 plotMultiStepDistances"
6. Long-term predictions for up to 100 steps for a case with increased inlet velocity u_in = 1.2 m/s:
   execute predictLongTerm.sh; the exact and the predicted mean flow field can be plotted with "python3 plotUFields.py"

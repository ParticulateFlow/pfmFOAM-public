IFS=$' ' read -a tlist < listOfReferenceStates

for t in "${tlist[@]}"
do
   printf "tstart %s;" "$t" > referenceState
   printf "mode internal;" > mode
   singlePhaseResponseFunctionCalculator > run_internal.log 2>&1
   printf "mode boundary;" > mode
   singlePhaseResponseFunctionCalculator > run_boundary.log 2>&1
   printf "mode integrated;" > mode
   singlePhaseResponseFunctionCalculator > run_integrated.log 2>&1
   mkdir dataBase/$t
   cp timeSeries/$t/U dataBase/$t/
   mv run*.log dataBase/$t/
   mv targetCells_* dataBase/$t/
   mv X_uu_* dataBase/$t/
done

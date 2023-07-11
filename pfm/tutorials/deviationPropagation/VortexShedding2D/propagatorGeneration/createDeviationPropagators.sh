IFS=$' ' read -a tlist < listOfReferenceStates
timeSeries=../referenceDataGeneration/data_referenceTimeSeries
mkdir dataBase

for t in "${tlist[@]}"
do
   printf "tstart %s;" "$t" > referenceState
   printf "mode internal;" > mode
   singlePhaseDeviationPropagatorCalculator > propagator_internal.log 2>&1
   printf "mode boundary;" > mode
   singlePhaseDeviationPropagatorCalculator > propagator_boundary.log 2>&1
   printf "mode integrated;" > mode
   singlePhaseDeviationPropagatorCalculator > propagator_integrated.log 2>&1
   mkdir dataBase/$t
   tend=$(<timeEvolvedRefState)
   cp $timeSeries/$t/U dataBase/$t/
   cp $timeSeries/${tend}/U dataBase/$t/UEvolved
   mv targetCells_* dataBase/$t/
   mv X_uu_* dataBase/$t/
done

#!/bin/bash
path="../validationData/uin_120/data_uin_120"
dt=1.0;
nsteps=20;


t3=0
for ((i=1; i<=$nsteps; i++))
do
# t3=$(awk '{print $1+$2}' <<<"${t3} ${dt}")
  t3=$(echo "$t3+$dt" | bc)
  # add a leading 0 for t3 < 1.0, remove trailing 0 and trailing .
  # in case of multiple trailing 0s, the second line would have to be repeated
  t3=$(echo "$t3" | sed 's/^\./0./')
  t3=$(echo "$t3" | sed 's/0$//')
  t3=$(echo "$t3" | sed 's/\.$//')
  mkdir $t3
done

mkdir log

for t in $path/*
do
  t=${t%*/}       # remove any trailing "/"
  t1="${t##*/}"   # save everything after the final "/"
  t2=${t1}
  solutions_exist=true
  for ((i=1; i<=$nsteps; i++))
  do
    echo $i
#   t2=$(awk '{print $1+$2}' <<<"${t2} ${dt}")
    t2=$(echo "$t2+$dt" | bc)
    t2=$(echo "$t2" | sed 's/^\./0./')
    t2=$(echo "$t2" | sed 's/0$//')
    t2=$(echo "$t2" | sed 's/\.$//')
    if [ -d "$path/$t2" ];
    then
      echo ""
    else
      solutions_exist=false
    fi
  done

  if [ "$solutions_exist" = true ] ;
  then
    rm -r 0
    mkdir 0
    echo "copying from $path/$t1 to 0"
    cp $path/$t1/U 0/
    # get p file to correct solution if necessary
    cp $path/$t1/p 0/
    rm distances
    t2=${t1}
    t3=0
    for ((i=1; i<=$nsteps; i++))
    do
#     t2=$(awk '{print $1+$2}' <<<"${t2} ${dt}")
#     t3=$(awk '{print $1+$2}' <<<"${t3} ${dt}")
      t2=$(echo "$t2+$dt" | bc)
      t2=$(echo "$t2" | sed 's/^\./0./')
      t2=$(echo "$t2" | sed 's/0$//')
      t2=$(echo "$t2" | sed 's/\.$//')

      t3=$(echo "$t3+$dt" | bc)
      t3=$(echo "$t3" | sed 's/^\./0./')
      t3=$(echo "$t3" | sed 's/0$//')
      t3=$(echo "$t3" | sed 's/\.$//')

      cd $t3 && rm *
      cd ..
      echo "copying from $path/$t2 to $t3"
      cp $path/$t2/U $t3/UExactSolution
    done

    # run simulation
    singlePhaseDeviationPropagationPredictor > log/run_$t1.log 2>&1
    cat distances >> distancesList

    mv distances distances_$t1
    echo "$path/$t1" >> testStateList
  fi
done

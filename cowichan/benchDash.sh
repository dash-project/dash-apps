#!/bin/bash

probIxStart=$1; probIxEnd=$2; jobsStart=$3; jobsEnd=$4; numberOfIterations=$5;
nRowsCols=40000

lP=(randmat thresh winnow outer product chain)
lI=(randmat_in thresh_in winnow_in outProd_in outProd_in)

if [[ $$CC == *icc ]] ; then module swap intel gnu &> /dev/null; fi;

echo $nRowsCols $nRowsCols 100 > thresh_in;\
echo $nRowsCols $nRowsCols $nRowsCols > winnow_in;\
echo $nRowsCols > outProd_in;\

for (( IX=$probIxStart ; IX<=$probIxEnd ; ++IX )); do
 
  for (( jobs=$jobsStart ; jobs<=$jobsEnd ; ++jobs )); do
    for (( it=$numberOfIterations ; it >= 1 ; --it )); do
      echo "run DASH ${lP[$IX]} with: $jobs jobs and $nRowsCols nRowsCols. Iterations left:$it"
      case $IX in
        0) echo $nRowsCols $nRowsCols $(( 1 + RANDOM % 666 )) | mpirun -n $jobs ./randmat/randmat --is_bench;;
        5) echo $nRowsCols $(( 1 + RANDOM % 666 )) 100 $nRowsCols | mpirun -n $jobs ./chain/chain --is_bench;;
        *) mpirun -n $jobs ./${lP[$IX]}/${lP[$IX]} ${lI[$IX]} --is_bench;;
      esac
    done
  done

done


#!/bin/bash
CHPL_ROOT=baSrcPaper/chapel

probIxStart=$1; probIxEnd=$2; jobsStart=$3; jobsEnd=$4; numberOfIterations=$5;
nRowsCols=4000

lP=(randmat thresh winnow outer product chain)
lI=(randmat_in thresh_in winnow_in outProd_in outProd_in)

if [[ $$CC == *icc ]] ; then module swap intel gnu &> /dev/null; fi;

echo $nRowsCols $nRowsCols 100 > thresh_in;\
echo $nRowsCols $nRowsCols $nRowsCols > winnow_in;\
echo $nRowsCols > outProd_in;\

for (( IX=$probIxStart ; IX<=$probIxEnd ; ++IX )); do
 
  for (( jobs=$jobsStart ; jobs<=$jobsEnd ; jobs+=4 )); do
    for (( it=$numberOfIterations ; it >= 1 ; --it )); do
      echo "run Chapel ${lP[$IX]} with: $jobs jobs and $nRowsCols nRowsCols. Iterations left:$it"
      case $IX in
        0) echo $nRowsCols $nRowsCols $(( 1 + RANDOM % 666 )) | $CHPL_ROOT/randmat/expertpar/main --dataParTasksPerLocale=$jobs --is_bench;;
        5) echo $nRowsCols $(( 1 + RANDOM % 666 )) 100 $nRowsCols | $CHPL_ROOT/chain/expertpar/main --dataParTasksPerLocale=$jobs --is_bench;;
        *) $CHPL_ROOT/${lP[$IX]}/expertpar/main --dataParTasksPerLocale=$jobs --is_bench < ${lI[$IX]};;
      esac
    done
  done

done


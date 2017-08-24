#!/bin/bash
TBB_ROOT=baSrcPaper/tbb

probIxStart=$1
probIxEnd=$2
jobsStart=$3
jobsEnd=$4
numberOfIterations=$5
nRowsCols=$6
thresh=$7
winnowNelts=$8

lP=(randmat thresh winnow outer product chain)
lI=(randmat_in thresh_in winnow_in outProd_in outProd_in)

if [[ $CC == *gcc ]] ; then module swap gnu intel &> /dev/null; fi;

echo $nRowsCols $nRowsCols $thresh > thresh_in
echo $nRowsCols $nRowsCols $winnowNelts $thresh > winnow_in
echo $nRowsCols > outProd_in

for (( IX=$probIxStart ; IX<=$probIxEnd ; ++IX )); do
 
  for (( jobs=$jobsStart ; jobs<=$jobsEnd ; jobs+=1 )); do
  
    for (( it=$numberOfIterations ; it >= 1 ; --it )); do
      echo "run TBB ${lP[$IX]} with: jobs:$jobs nRowsCols:$nRowsCols thresh:$thresh win_nelts:$winnowNelts ItsLeft:$it"
      case $IX in
        0) echo $nRowsCols $nRowsCols $(( 1 + RANDOM % 666 )) | $TBB_ROOT/randmat/expertpar/main --threads $jobs --is_bench;;
        5) echo $nRowsCols $(( 1 + RANDOM % 666 )) $thresh $winnowNelts | $TBB_ROOT/chain/expertpar/main --threads $jobs --is_bench;;
        *) $TBB_ROOT/${lP[$IX]}/expertpar/main --threads $jobs --is_bench < ${lI[$IX]};;
      esac
    done
  done

done


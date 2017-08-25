#!/bin/bash
GO_ROOT=baSrcPaper/go

probIxStart=$1
probIxEnd=$2
jobsStart=$3
jobsEnd=$4
numberOfIterations=$5
nRowsCols=$6
thresh=$7
winnowNelts=$8

CYAN='\033[0;36m'
GREEN='\033[0;32m'
NC='\033[0m'
BCOLORS=('\033[0;41m' '\033[0;42m' '\033[0;43m' '\033[0;44m' '\033[0;45m' '\033[0;46m')

lP=(randmat thresh winnow outer product chain)
lI=(randmat_in thresh_in winnow_in outProd_in outProd_in)

if [[ $CC == *icc ]] ; then module swap intel gnu &> /dev/null; fi;

echo $nRowsCols $nRowsCols $thresh > thresh_in
echo $nRowsCols $nRowsCols $winnowNelts $thresh > winnow_in
echo $nRowsCols > outProd_in

for (( IX=$probIxStart ; IX<=$probIxEnd ; ++IX )); do
 
  for (( jobs=$jobsStart ; jobs<=$jobsEnd ; jobs+=1 )); do
    export GOMAXPROCS=$jobs
    for (( it=$numberOfIterations ; it >= 1 ; --it )); do
      # echo "run Go ${lP[$IX]} with: jobs:$jobs nRowsCols:$nRowsCols thresh:$thresh win_nelts:$winnowNelts ItsLeft:$it"
      printf "run...$GREEN GO$NC     ${BCOLORS[$IX]}%7s$NC->jobs:$CYAN$jobs$NC nRowsCols:$CYAN%5u$NC thresh:$CYAN%3u$NC win_nelts:$CYAN%5u$NC ItsLeft:$CYAN$it$NC\n" ${lP[$IX]} $nRowsCols $thresh $winnowNelts
      case $IX in
        0) echo $nRowsCols $nRowsCols $(( 1 + RANDOM % 666 )) | $GO_ROOT/randmat/expertpar/main --is_bench;;
        5) echo $nRowsCols $(( 1 + RANDOM % 666 )) $thresh $winnowNelts | $GO_ROOT/chain/expertpar/main --is_bench;;
        *) $GO_ROOT/${lP[$IX]}/expertpar/main --is_bench < ${lI[$IX]};;
      esac
    done
  done

done


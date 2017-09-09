#!/bin/bash
CILK_ROOT=baSrcPaper/cilk

probIxStart=$1
probIxEnd=$2
jobsStart=$3
jobsEnd=$4
numberOfIterations=$5
nRowsCols=$6
thresh=$7
winnowNelts=$8

CYAN='\033[0;36m'
BLUE='\033[0;34m'
NC='\033[0m'
BCOLORS=('\033[0;41m' '\033[0;42m' '\033[0;43m' '\033[0;44m' '\033[0;45m' '\033[0;46m')

lP=(randmat thresh winnow outer product chain)
lI=(randmat_in thresh_in winnow_in outProd_in outProd_in)

if [[ $CC == *gcc ]] ; then module swap gnu intel &> /dev/null; fi;

echo $nRowsCols $nRowsCols $thresh > thresh_in
echo $nRowsCols $nRowsCols $winnowNelts $thresh > winnow_in
echo $winnowNelts > outProd_in

for (( IX=$probIxStart ; IX<=$probIxEnd ; ++IX )); do
 
  for (( jobs=$jobsStart ; jobs<=$jobsEnd ; jobs+=1 )); do
    export CILK_NWORKERS=$jobs
    for (( it=$numberOfIterations ; it >= 1 ; --it )); do
      # echo -e "run$BLUE Cilk$NC ${BCOLORS[$IX]}${lP[$IX]}$NC with: jobs:$CYAN$jobs$NC nRowsCols:$CYAN$nRowsCols$NC thresh:$CYAN$thresh$NC win_nelts:$CYAN$winnowNelts$NC ItsLeft:$CYAN$it$NC"
      printf "run...$BLUE Cilk$NC   ${BCOLORS[$IX]}%7s$NC->jobs:$CYAN%2u$NC nRowsCols:$CYAN%5u$NC thresh:$CYAN%3u$NC win_nelts:$CYAN%5u$NC ItsLeft:$CYAN$it$NC\n" ${lP[$IX]} $jobs $nRowsCols $thresh $winnowNelts
      case $IX in
        0) echo $nRowsCols $nRowsCols $(( 1 + RANDOM % 666 )) | $CILK_ROOT/randmat/expertpar/main --nproc $jobs --is_bench;;
        5) echo $nRowsCols $(( 1 + RANDOM % 666 )) $thresh $winnowNelts | $CILK_ROOT/chain/expertpar/main --nproc $jobs --is_bench;;
        *) $CILK_ROOT/${lP[$IX]}/expertpar/main --nproc $jobs --is_bench < ${lI[$IX]};;
      esac
    done
  done

done


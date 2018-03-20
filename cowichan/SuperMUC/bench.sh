#!/bin/bash

CYAN='\033[0;36m'
RED='\033[0;31m'
NC='\033[0m'
BCOLORS=('\033[0;41m' '\033[0;42m' '\033[0;43m' '\033[0;44m' '\033[0;45m' '\033[0;46m')

NodesStart=1
NodesEnd=1
Procs=28
Threads=1
numberOfRuns=1
nRowsCols=(100 400 4000 40000)
thresh=(50 100)
winnowNelts=(100 400 4000 40000)


for (( Nodes=$NodesStart ; Nodes<=$NodesEnd ; ++Nodes )); do
  for (( run=1 ; run <= $numberOfRuns ; ++run )); do



    # #bench randmat
    for (( nRowsColsIX=1 ; nRowsColsIX<4 ; ++nRowsColsIX )); do
      printf "submit...\033[0;41m%7s$NC->nodes$CYAN%2u$NC procs$CYAN%2u$NC threads$CYAN%2u$NC mx:$CYAN%6u$NC run:$CYAN$run$NC\n" \
              randmat $Nodes $Procs $Threads ${nRowsCols[$nRowsColsIX]}
      benchrandmat ~/template ${nRowsCols[$nRowsColsIX]} $run $Nodes $Threads $Procs
      break
    done

    # #bench thresh
    for (( nRowsColsIX=1 ; nRowsColsIX<4 ; ++nRowsColsIX )); do
      for (( threshIX=1 ; threshIX<2 ; ++threshIX )); do
        printf "submit...\033[0;42m%7s$NC->nodes$CYAN%2u$NC procs$CYAN%2u$NC threads$CYAN%2u$NC mx:$CYAN%6u$NC trsh:$CYAN%3u$NC run:$CYAN$run$NC\n" \
                thresh $Nodes $Procs $Threads ${nRowsCols[$nRowsColsIX]} ${thresh[$threshIX]}
        benchthresh ~/template ${nRowsCols[$nRowsColsIX]} ${thresh[$threshIX]} $run $Nodes $Threads $Procs
        break 2
      done
    done

    #bench winnows and chains
    for (( nRowsColsIX=1 ; nRowsColsIX<4 ; ++nRowsColsIX )); do
      for (( threshIX=1 ; threshIX<2 ; ++threshIX )); do
        for (( winnowNeltsIX=1 ; winnowNeltsIX<4 ; ++winnowNeltsIX )); do
          #bench winnow
          printf "submit...\033[0;43m%7s$NC->nodes$CYAN%2u$NC procs$CYAN%2u$NC threads$CYAN%2u$NC mx:$CYAN%6u$NC trsh:$CYAN%3u$NC nelts:$CYAN%3u$NC run:$CYAN$run$NC\n" \
                winnow $Nodes $Procs $Threads ${nRowsCols[$nRowsColsIX]} ${thresh[$threshIX]} ${winnowNelts[$winnowNeltsIX]}
          benchwinnow ~/template ${nRowsCols[$nRowsColsIX]} ${thresh[$threshIX]} ${winnowNelts[$winnowNeltsIX]} $run $Nodes $Threads $Procs

          # bench chains
          printf "submit...\033[0;46m%7s$NC->nodes$CYAN%2u$NC procs$CYAN%2u$NC threads$CYAN%2u$NC mx:$CYAN%6u$NC trsh:$CYAN%3u$NC nelts:$CYAN%3u$NC run:$CYAN$run$NC\n" \
                chain $Nodes $Procs $Threads ${nRowsCols[$nRowsColsIX]} ${thresh[$threshIX]} ${winnowNelts[$winnowNeltsIX]}
          benchchain ~/template ${nRowsCols[$nRowsColsIX]} ${thresh[$threshIX]} ${winnowNelts[$winnowNeltsIX]} $run $Nodes $Threads $Procs
          break 3
        done
      done
    done

    #bench outers and products
    for (( winnowNeltsIX=1 ; winnowNeltsIX<4 ; ++winnowNeltsIX )); do
        printf "submit...\033[0;44m%7s$NC->nodes$CYAN%2u$NC procs$CYAN%2u$NC threads$CYAN%2u$NC nelts:$CYAN%6u$NC run:$CYAN$run$NC\n" \
              outer $Nodes $Procs $Threads ${winnowNelts[$winnowNeltsIX]}
        benchouter ~/template ${winnowNelts[$winnowNeltsIX]} $run $Nodes $Threads $Procs

        printf "submit...\033[0;45m%7s$NC->nodes$CYAN%2u$NC procs$CYAN%2u$NC threads$CYAN%2u$NC nelts:$CYAN%6u$NC run:$CYAN$run$NC\n" \
              product $Nodes $Procs $Threads ${winnowNelts[$winnowNeltsIX]}
        benchproduct ~/template ${winnowNelts[$winnowNeltsIX]} $run $Nodes $Threads $Procs
        break
    done


  done
done


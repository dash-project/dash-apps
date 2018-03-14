#!/bin/bash


probIxStart=0
probIxEnd=5
jobsStart=1
jobsEnd=24
numberOfIterations=10
nRowsCols=(100 400 4000 40000)
thresh=(10 25 50 75 100)
winnowNelts=(100 400 4000 40000)


# #bench winnows and chains
# for (( nRowsColsIX=3 ; nRowsColsIX<4 ; ++nRowsColsIX )); do
  # for (( threshIX=3 ; threshIX<5 ; ++threshIX )); do
    # for (( winnowNeltsIX=0 ; winnowNeltsIX<4 ; ++winnowNeltsIX )); do
      # #bench winnows
      # ./benchDash.sh 2 2 $jobsStart $jobsEnd $numberOfIterations ${nRowsCols[$nRowsColsIX]} ${thresh[$threshIX]} ${winnowNelts[$winnowNeltsIX]}
      # #bench chains
      # ./benchDash.sh 5 5 $jobsStart $jobsEnd $numberOfIterations ${nRowsCols[$nRowsColsIX]} ${thresh[$threshIX]} ${winnowNelts[$winnowNeltsIX]}
    # done
  # done
# done


# #bench randmats
# for (( nRowsColsIX=2 ; nRowsColsIX<4 ; ++nRowsColsIX )); do
  # ./benchDash.sh 0 0 $jobsStart $jobsEnd $numberOfIterations ${nRowsCols[$nRowsColsIX]} 0 0
# done

# #bench threshs
# tmp=2
# for (( nRowsColsIX=2 ; nRowsColsIX<4 ; ++nRowsColsIX )); do
  # for (( threshIX=tmp ; threshIX<5 ; ++threshIX )); do
  # tmp=0
    # ./benchDash.sh 1 1 $jobsStart $jobsEnd $numberOfIterations ${nRowsCols[$nRowsColsIX]} ${thresh[$threshIX]} 0
  # done
# done


#bench outers and products
for (( winnowNeltsIX=0 ; winnowNeltsIX<4 ; ++winnowNeltsIX )); do
    ./benchDash.sh 3 4 $jobsStart $jobsEnd $numberOfIterations 0 0 ${winnowNelts[$winnowNeltsIX]}
done
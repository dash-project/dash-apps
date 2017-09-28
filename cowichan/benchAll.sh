#!/bin/bash


probIxStart=0
probIxEnd=5
jobsStart=1
jobsEnd=24
numberOfIterations=5
nRowsCols=(100 400 4000 40000)
thresh=(10 25 50 75 100)
winnowNelts=(100 400 4000 40000)


for (( nRowsColsIX=0 ; nRowsColsIX<4 ; ++nRowsColsIX )); do
  for (( threshIX=0 ; threshIX<5 ; ++threshIX )); do
    for (( winnowNeltsIX=0 ; winnowNeltsIX<4 ; ++winnowNeltsIX )); do
      ./benchDash.sh $probIxStart $probIxEnd $jobsStart $jobsEnd $numberOfIterations ${nRowsCols[$nRowsColsIX]} ${thresh[$threshIX]} ${winnowNelts[$winnowNeltsIX]}
      ./benchTBB.sh $probIxStart $probIxEnd $jobsStart $jobsEnd $numberOfIterations ${nRowsCols[$nRowsColsIX]} ${thresh[$threshIX]} ${winnowNelts[$winnowNeltsIX]}
      ./benchChapel.sh $probIxStart $probIxEnd $jobsStart $jobsEnd $numberOfIterations ${nRowsCols[$nRowsColsIX]} ${thresh[$threshIX]} ${winnowNelts[$winnowNeltsIX]}
      ./benchCilk.sh $probIxStart $probIxEnd $jobsStart $jobsEnd $numberOfIterations ${nRowsCols[$nRowsColsIX]} ${thresh[$threshIX]} ${winnowNelts[$winnowNeltsIX]}
      ./benchGo.sh $probIxStart $probIxEnd $jobsStart $jobsEnd $numberOfIterations ${nRowsCols[$nRowsColsIX]} ${thresh[$threshIX]} ${winnowNelts[$winnowNeltsIX]}
    done
  done
done
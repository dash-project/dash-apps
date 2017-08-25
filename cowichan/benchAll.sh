#!/bin/bash


probIxStart=0
probIxEnd=5
jobsStart=24
jobsEnd=24
numberOfIterations=1
nRowsCols=(100 400 4000 40000)
thresh=(10 25 50 75 100)
winnowNelts=(100 400 4000 40000)


for (( nRowsColsIX=0 ; nRowsColsIX<1 ; ++nRowsColsIX )); do
  for (( threshIX=0 ; threshIX<1 ; ++threshIX )); do
    for (( winnowNeltsIX=0 ; winnowNeltsIX<1 ; ++winnowNeltsIX )); do
      ./benchDash.sh $probIxStart $probIxEnd $jobsStart $jobsEnd $numberOfIterations ${nRowsCols[$nRowsColsIX]} ${thresh[$threshIX]} ${winnowNelts[$winnowNeltsIX]}
      ./benchTBB.sh $probIxStart $probIxEnd $jobsStart $jobsEnd $numberOfIterations ${nRowsCols[$nRowsColsIX]} ${thresh[$threshIX]} ${winnowNelts[$winnowNeltsIX]}
      ./benchChapel.sh $probIxStart $probIxEnd $jobsStart $jobsEnd $numberOfIterations ${nRowsCols[$nRowsColsIX]} ${thresh[$threshIX]} ${winnowNelts[$winnowNeltsIX]}
      ./benchCilk.sh $probIxStart $probIxEnd $jobsStart $jobsEnd $numberOfIterations ${nRowsCols[$nRowsColsIX]} ${thresh[$threshIX]} ${winnowNelts[$winnowNeltsIX]}
      ./benchGo.sh $probIxStart $probIxEnd $jobsStart $jobsEnd $numberOfIterations ${nRowsCols[$nRowsColsIX]} ${thresh[$threshIX]} ${winnowNelts[$winnowNeltsIX]}
    done
  done
done
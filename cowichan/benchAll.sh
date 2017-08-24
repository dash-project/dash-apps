#!/bin/bash


probIxStart=0
probIxEnd=5
jobsStart=4
jobsEnd=4
numberOfIterations=1
nRowsCols=20
thresh=50
winnowNelts=20


./benchDash.sh $probIxStart $probIxEnd $jobsStart $jobsEnd $numberOfIterations $nRowsCols $thresh $winnowNelts
./benchTBB.sh $probIxStart $probIxEnd $jobsStart $jobsEnd $numberOfIterations $nRowsCols $thresh $winnowNelts
./benchChapel.sh $probIxStart $probIxEnd $jobsStart $jobsEnd $numberOfIterations $nRowsCols $thresh $winnowNelts
./benchCilk.sh $probIxStart $probIxEnd $jobsStart $jobsEnd $numberOfIterations $nRowsCols $thresh $winnowNelts
./benchGo.sh $probIxStart $probIxEnd $jobsStart $jobsEnd $numberOfIterations $nRowsCols $thresh $winnowNelts

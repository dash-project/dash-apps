#!/bin/bash

for (( jobs=1; jobs <= 24 ; ++jobs )); do
  echo "lauf mit: $jobs jobs"
  mpirun -n $jobs ./winnow/winnow "winnow_in" --is_bench
done
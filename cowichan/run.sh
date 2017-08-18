#!/bin/bash

for (( jobs=1; jobs <= 24 ; ++jobs )); do
  echo "run with: $jobs jobs"
  mpirun -n $jobs ./outer/outer "outProd_in" --is_bench
done
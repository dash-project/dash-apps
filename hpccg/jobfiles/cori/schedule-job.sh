#!/bin/bash

[[ !  -z  "$1"  ]] || { echo "usage: $0 <Jobtemplate>"; exit 1; }
TPL="$1"
UNITS="1 2 4 8 16 32"

for NUNITS in $UNITS; do
  JOBOUT="$TPL.$NUNITS"
  sed "s/<procs>/$NUNITS/g" $TPL > $JOBOUT 
  sbatch $JOBOUT
  rm $JOBOUT 
done


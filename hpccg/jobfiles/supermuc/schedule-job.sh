#!/bin/bash

[[ !  -z  "$1"  ]] || { echo "usage: $0 <Jobtemplate>"; exit 1; }
TPL="$1"
UNITS="1 2 4 8 16 28 56 112 224"
PPN="28"

for NUNITS in $UNITS; do
  JOBOUT="$TPL.$NUNITS"
  NODES=$((($NUNITS-1) / ($PPN+1) +1))
  sed "s/<procs>/$NUNITS/g" $TPL > $JOBOUT
  sed -i "s/<nodes>/$NODES/g" $JOBOUT
  llsubmit $JOBOUT
  rm $JOBOUT
done


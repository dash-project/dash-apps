#!/bin/bash

[[ !  -z  "$1"  ]] || { echo "usage: $0 <ENVIRONMENT>  for building Makefile.ENVIRONMENT"; exit 1; }
ENVIRONMENT="$1"

for DIR in $(ls -d */); do
  if [[ "$DIR" == "jobfiles/" ]];then continue; fi
  if [[ "$DIR" == "results/" ]];then continue; fi
  test -e $DIR/Makefile.$ENVIRONMENT || { echo "no Makefile for $ENVIRONMENT in $DIR"; continue; }

  cd $DIR
  echo "Build $DIR"
  make clean
  make -f Makefile.$ENVIRONMENT -j 10
  cd ..
done

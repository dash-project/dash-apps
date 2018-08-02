#!/bin/bash

if [[ ! -f "$1" || ! -f "$2" ]]
then
  echo "usage: $0 <left> <right>"
  exit 1
fi

colordiff --strip-trailing-cr <(sed -n 19,125p "$1" | sed -r 's/\s+//g') <(sed -n 19,125p "$2" | sed -r 's/\s+//g')


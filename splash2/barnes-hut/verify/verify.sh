#!/bin/bash

if [[ ! -f "$1" || ! -f "$2" ]]
then
  echo -e "usage: $0 <lhs> <rhs>\n\
\
NOTE: This requires the output of the final tree. See 'printtree' in the source code."
  exit 1
fi

gpat='^Cell\|Child\|Leaf\|Body'
spat='s/ //g;s/Num=[0-9]\+,//g'

diff -w \
  <(grep "$gpat" "$1" | sed "$spat" | sort) \
  <(grep "$gpat" "$2" | sed "$spat" | sort)


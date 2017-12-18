#!/bin/bash

# if present, the first command line argument is used as a label. It should not contain spaces.


cat trace_0*.csv >"trace"$1".csv"
rm -Rf trace_0*.csv

shopt -s nullglob
for i in image_unit0.csv.*; do
    NUM=`echo $i|sed -e 's/image_unit0.csv.//'`
    echo -n $NUM" "
    #cat "image_unit"*".csv."$NUM | sort > "image.csv."$NUM
    sort -n -m "image_unit"*".csv."$NUM > "image_"$1".csv."$NUM
    rm "image_unit"*".csv."$NUM
done

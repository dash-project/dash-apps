#!/bin/bash

for i in image_unit0.csv.*; do
    NUM=`echo $i|sed -e 's/image_unit0.csv.//'`
    echo -n $NUM" "
    #cat "image_unit"*".csv."$NUM | sort > "image.csv."$NUM
    sort -m "image_unit"*".csv."$NUM > "image.csv."$NUM
    rm "image_unit"*".csv."$NUM
done

cat trace0*.csv >trace_.csv
rm -Rf trace0*.csv

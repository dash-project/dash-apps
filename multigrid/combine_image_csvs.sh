#!/bin/bash

for i in image_unit0.csv.*; do
    NUM=`echo $i|sed -e 's/image_unit0.csv.//'`
    echo -n $NUM" "
    cat "image_unit"*".csv."$NUM | sort > "image.csv."$NUM
    rm "image_unit"*".csv."$NUM
done

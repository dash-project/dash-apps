#!/bin/bash

# if present, the first command line argument is used as a label. It should not contain spaces.

cat overview_0*.csv >"overview"$1".csv"
rm -Rf overview_0*.csv
`dirname $0`/overview.gnuplot

ls image_unit0.csv.* 2>/dev/null || exit 0

for i in image_unit0.csv.*; do
    NUM=`echo $i|sed -e 's/image_unit0.csv.//'`
    echo -n $NUM" "
    ##echo sort -n -m "image_unit"*".csv."$NUM ">" "image_"$1".csv."$NUM
    sort -n -m "image_unit"*".csv."$NUM > "image_"$1".csv."$NUM
    rm "image_unit"*".csv."$NUM
done


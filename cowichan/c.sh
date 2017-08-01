#!/bin/bash

# findet alle mains unter outer/expertpar
# out=($(find . -regextype posix-extended -regex "\./.*/outer/expertpar/main\.(go|[chpl]{2,4})"))
# echo ${out[2]}

#how to use:
# cow && cd $( ./c.sh thresh 0)
# -> leads to folder: ./baSrcPaper/chapel/thresh/expertpar
l=(chapel go cilk tbb)

echo ./baSrcPaper/${l[$2]}/$1/expertpar


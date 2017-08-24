#!/bin/bash

# findet alle mains unter outer/expertpar
# out=($(find . -regextype posix-extended -regex "\./.*/outer/expertpar/main\.(go|[chpl]{2,4})"))
# echo ${out[2]}

#how to use:
# cow && cd $( ./c.sh thresh 0)
# -> leads to folder: ./baSrcPaper/chapel/thresh/expertpar
lL=(chapel go cilk tbb)
lP=(randmat thresh winnow outer product chain)

echo ./baSrcPaper/${lL[$1]}/${lP[$2]}/expertpar


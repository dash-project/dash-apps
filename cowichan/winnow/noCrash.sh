#!/bin/bash

make && mpirun -n 1 ./winnow "winnow_in3" --is_bench
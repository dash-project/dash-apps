#!/bin/bash

make && mpirun -n 1 ./winnow "winnow_in2" --is_bench
#!/bin/bash

make && mpirun -n 2 ./winnow "winnow_in2" --is_bench
#!/bin/bash

make && mpirun -n 2 ./winnow "winnow_in3" --is_bench
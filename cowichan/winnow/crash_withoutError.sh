#!/bin/bash

make && mpirun -n 1 ./winnow "winnow_in" --is_bench
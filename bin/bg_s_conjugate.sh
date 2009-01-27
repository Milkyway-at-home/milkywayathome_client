#!/bin/bash

echo 'job starting'
mpirun -mode VN -cwd `pwd` ./bg_astronomy -synch -stars stars-86.txt -parameters parameters-86-modified.txt -cgd -gd_reset 8 -gd_min_threshold 0.000001 -gd_iterations 5
# additional calls to mpirun may be placed in the script, they will all use the same partition

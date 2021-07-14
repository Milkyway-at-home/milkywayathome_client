#!/bin/sh
cd ..
make gpu_test -j16
cd bin/
./gpu_test
python3 read.py "test_outputCPU.txt" "test_outputCPU.dat"
python3 read.py "test_outputGPU.txt" "test_outputGPU.dat"
# python3 read.py "output0PU.txt" "test_output0PU.dat"
# python3 read.py "output0PU.txt" "test_output0PU.dat"
gnuplot << EOF
set terminal png
set output "trial.png"
plot "test_outputCPU.dat", "test_outputGPU.dat"
EOF
python3 error.py "test_outputGPU.dat" "test_outputCPU.dat" "error.dat"

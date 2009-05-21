nvcc -lcuda  -I/Developer/CUDA/common/inc -L/Developer/CUDA/common/lib likelihood.cu *.c ../util/io_util.c ../astronomy/parameters.c ../astronomy/star_points.c -o likelihood_test

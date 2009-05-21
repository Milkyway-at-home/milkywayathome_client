g++ -I/usr/local/cuda/include -c -o coords.o coords.c
g++ -I/usr/local/cuda/include -c -o gauss_legendre.o gauss_legendre.c
g++ -I/usr/local/cuda/include -c -o cpu_r_constants.o cpu_r_constants.c
g++ -I/usr/local/cuda/include -c -o cpu_coords.o cpu_coords.c
g++ -I/usr/local/cuda/include -c -o cpu_integrals.o cpu_integrals.c
nvcc -lcuda -I/Developer/CUDA/common/inc -L/usr/local/cuda/lib -c -o evaluation_gpu.o evaluation_gpu.cu
g++ -DBOINC_APPLICATION -DBOINC_APP_NAME='"stock_cuda_gpu"' -DBOINC_APP_VERSION=0.01 -DMILKYWAY_GPU -DCOMPUTE_ON_GPU \
	-Wall \
	-I/Developer/CUDA/common/inc -I/usr/local/cuda/include -L/usr/local/cuda/lib -lcudart \
	-I/Users/deselt/software/boinc -I/Users/deselt/software/boinc/clientgui/mac -I/Users/deselt/software/boinc/api -I/Users/deselt/software/boinc/lib -L/Users/deselt/software/boinc/mac_build/build/Deployment -lboinc -lboinc_api \
	../astronomy/boinc_astronomy.c ../astronomy/parameters.c ../astronomy/star_points.c ../astronomy/evaluation_state.c \
	../util/io_util.c ../util/matrix.c \
	../searches/hessian.c ../searches/gradient.c ../searches/line_search.c ../searches/newton_method.c ../searches/regression.c ../searches/search_parameters.c ../searches/result.c \
	../evaluation/simple_evaluator.c \
	gauss_legendre.o coords.o cpu_r_constants.o cpu_coords.o cpu_integrals.o \
	evaluation_gpu.o \
	-o milkyway_gpu_0.01_x86_64-apple-darwin

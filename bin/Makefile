APP_VERSION = 0.19
APP_NAME = your_app_name-gpu-test

OSX_BOINC_DIR = /Users/deselt/software/boinc
OSX_LIBBOINC_DIR = /Users/deselt/software/boinc/mac_build/build/Deployment

LINUX_BOINC_DIR = /software/boinc-6.3.14-0/pkg/include/boinc
LINUX_LIBBOINC_DIR = /usr/local/lib
#/software/boinc-6.3.14-0/pkg/lib

CUDA_INC = -I/Developer/CUDA/common/inc -I/usr/local/cuda/include -I /opt/cuda/sdk/common/inc
CUDA_LIB = -L/usr/local/cuda/lib -L/opt/cuda/sdk/lib -L/opt/cuda/lib


GPU_FLAGS = -DMILKYWAY_GPU -DCOMPUTE_ON_GPU $(CUDA_INC) $(CUDA_LIB)
#GPU_PRECISION = -DDOUBLE_PRECISION
GPU_PRECISION = -DSINGLE_PRECISION
CPU_FLAGS = -DMILKYWAY -DCOMPUTE_ON_CPU

VARIANTFLAGS = -DBOINC_APPLICATION -DBOINC_APP_NAME='"$(APP_NAME)"' -DBOINC_APP_VERSION=$(APP_VERSION)
CXXFLAGS_ALL = $(VARIANTFLAGS) -g


OSX_CC_PPC = /usr/bin/gcc-3.3
OSX_CXX_PPC = /usr/bin/g++-3.3
OSX_CXXFLAGS_PPC = -arch ppc -D_NONSTD_SOURCE -isystem /Developer/SDKs/MacOSX10.3.9.sdk -O2 -maltivec -mabi=altivec -mcpu=7400 -funroll-loops -DMAC_OS_X_VERSION_MIN_REQUIRED=1030 $(CXXFLAGS_ALL) \
	-I$(OSX_BOINC_DIR) -I$(OSX_BOINC_DIR)/api -I$(OSX_BOINC_DIR)/lib -I$(OSX_BOINC_DIR)/clientgui/mac/ \
	-L$(OSX_BOINC_DIR)/mac_build/build/Deployment
OSX_LDFLAGS_PPC = -Wl,-syslibroot,/Developer/SDKs/MacOSX10.3.9.sdk,-arch,ppc

OSX_CC_i686 = /usr/bin/gcc-4.0
OSX_CXX_i686 = /usr/bin/g++-4.0
OSX_CXXFLAGS_i686 = -arch i386 -isysroot /Developer/SDKs/MacOSX10.4u.sdk -mmacosx-version-min=10.4 -O3 -msse3 -mfpmath=sse -mtune=prescott -ftree-vectorize -funroll-loops $(CXXFLAGS_ALL) \
	-I$(OSX_BOINC_DIR) -I$(OSX_BOINC_DIR)/api -I$(OSX_BOINC_DIR)/lib -I$(OSX_BOINC_DIR)/clientgui/mac/ \
	-L$(OSX_BOINC_DIR)/mac_build/build/Deployment
OSX_LDFLAGS_i686 = -Wl,-syslibroot,/Developer/SDKs/MacOSX10.4u.sdk,-arch,i386

OSX_CC_x86_64 = /usr/bin/gcc-4.0
OSX_CXX_x86_64 = /usr/bin/g++-4.0
OSX_CXXFLAGS_x86_64 = -arch x86_64 -isysroot /Developer/SDKs/MacOSX10.5.sdk -mmacosx-version-min=10.5 -O3 -mfpmath=sse -march=nocona -msse3 -ftree-vectorize -funroll-loops $(CXXFLAGS_ALL) \
	-I$(OSX_BOINC_DIR) -I$(OSX_BOINC_DIR)/api -I$(OSX_BOINC_DIR)/lib -I$(OSX_BOINC_DIR)/clientgui/mac/ \
	-L$(OSX_BOINC_DIR)/mac_build/build/Deployment
OSX_LDFLAGS_x86_64 = -Wl,-syslibroot,/Developer/SDKs/MacOSX10.5.sdk,-arch,x86_64


LINUX_CXX_i686 = g++
LINUX_CXXFLAGS_i686 = -m32 -O2 -ftree-vectorize -funroll-loops $(CXXFLAGS_ALL) \
	-I$(LINUX_BOINC_DIR) -L$(LINUX_LIBBOINC_DIR)
LINUX_LDFLAGS_i686 = -L/usr/X11R6/lib -L/usr/local/lib

LINUX_CXX_x86_64 = g++
LINUX_CXXFLAGS_x86_64 = -O2 -msse2 -ftree-vectorize -funroll-loops $(CXXFLAGS_ALL) \
	-I$(LINUX_BOINC_DIR) -L$(LINUX_LIBBOINC_DIR) -I/usr/include/boinc
LINUX_LDFLAGS_x86_64 = -L/usr/local/lib -L/usr/lib64


FGDO_DIR = ..
APP_DIR = $(FGDO_DIR)/astronomy
GPU_DIR = $(FGDO_DIR)/astronomy_gpu
EVAL_DIR = $(FGDO_DIR)/evaluation
SEARCH_DIR = $(FGDO_DIR)/searches
UTIL_DIR = $(FGDO_DIR)/util

APP_OBJS = \
	$(APP_DIR)/boinc_astronomy.o \
	$(APP_DIR)/atSurveyGeometry.o \
	$(APP_DIR)/numericalIntegration.o \
	$(APP_DIR)/parameters.o \
	$(APP_DIR)/probability.o \
	$(APP_DIR)/stCoords.o \
	$(APP_DIR)/stCnum.o \
	$(APP_DIR)/stMath.o \
	$(APP_DIR)/stVector.o \
	$(APP_DIR)/star_points.o \
	$(APP_DIR)/evaluation_optimized.o \
	$(APP_DIR)/evaluation_state.o
SEARCH_OBJS = \
	$(SEARCH_DIR)/search_parameters.o \
	$(SEARCH_DIR)/result.o
UTIL_OBJS = \
	$(UTIL_DIR)/matrix.o \
	$(UTIL_DIR)/io_util.o \
	$(UTIL_DIR)/settings.o
GPU_APP_OBJS = \
	$(GPU_DIR)/evaluation_gpu6.o \
	$(GPU_DIR)/coords.o \
	$(GPU_DIR)/gauss_legendre.o \
	$(GPU_DIR)/cpu_r_constants.o \
	$(GPU_DIR)/cpu_coords.o \
	$(GPU_DIR)/cpu_integrals.o \
	$(EVAL_DIR)/simple_evaluator.o \
	$(SEARCH_DIR)/hessian.o \
	$(SEARCH_DIR)/gradient.o \
	$(SEARCH_DIR)/newton_method.o \
	$(SEARCH_DIR)/line_search.o \
	$(SEARCH_DIR)/regression.o

PROGS = milkyway_$(APP_VERSION)_powerpc-apple-darwin milkyway_$(APP_VERSION)_i686-apple-darwin milkyway_$(APP_VERSION)_x86_64-apple-darwin

all: $(PROGS)

$(GPU_DIR)/evaluation_gpu.o:
	nvcc -lcuda $(CUDA_INC) $(CUDA_LIB) -c -o $(GPU_DIR)/evaluation_gpu.o $(GPU_DIR)/evaluation_gpu.cu
$(GPU_DIR)/evaluation_gpu2.o: $(GPU_DIR)/evaluation_gpu2.cu
	nvcc -lcuda $(CUDA_INC) $(CUDA_LIB) -c -o $(GPU_DIR)/evaluation_gpu2.o $(GPU_DIR)/evaluation_gpu2.cu
$(GPU_DIR)/evaluation_gpu3.o:
	nvcc -lcuda $(GPU_PRECISION) $(CUDA_INC) $(CUDA_LIB) -c -o $(GPU_DIR)/evaluation_gpu3.o $(GPU_DIR)/evaluation_gpu3.cu
$(GPU_DIR)/evaluation_gpu4.o:
	nvcc -lcuda $(GPU_PRECISION) $(CUDA_INC) $(CUDA_LIB) -c -o $(GPU_DIR)/evaluation_gpu4.o $(GPU_DIR)/evaluation_gpu4.cu
$(GPU_DIR)/evaluation_gpu5.o:
	nvcc -lcuda $(GPU_PRECISION) $(CUDA_INC) $(CUDA_LIB) -c -o $(GPU_DIR)/evaluation_gpu5.o $(GPU_DIR)/evaluation_gpu5.cu
$(GPU_DIR)/evaluation_gpu6.o:
	nvcc -lcuda $(GPU_PRECISION) $(CUDA_INC) $(CUDA_LIB) -c -o $(GPU_DIR)/evaluation_gpu6.o $(GPU_DIR)/evaluation_gpu6.cu

osx_ppc: OBJ_CXX = $(OSX_CXX_PPC)
osx_ppc: OBJ_CXXFLAGS = $(CPU_FLAGS) $(OSX_CXXFLAGS_PPC)
osx_ppc: $(APP_OBJS) $(SEARCH_OBJS) $(UTIL_OBJS) $(OSX_LIBBOINC_DIR)/libboinc_api.a $(OSX_LIBBOINC_DIR)/libboinc.a
	$(OBJ_CXX) $(OBJ_CXXFLAGS) $(OSX_LDFLAGS_PPC) -o milkyway_$(APP_VERSION)_powerpc-apple-darwin $(APP_OBJS) $(SEARCH_OBJS) $(UTIL_OBJS) -lboinc -lboinc_api

osx_i686: OBJ_CXX = $(OSX_CXX_i686)
osx_i686: OBJ_CXXFLAGS = $(CPU_FLAGS) $(OSX_CXXFLAGS_i686)
osx_i686: $(APP_OBJS) $(SEARCH_OBJS) $(UTIL_OBJS) $(OSX_LIBBOINC_DIR)/libboinc_api.a $(OSX_LIBBOINC_DIR)/libboinc.a
	$(OBJ_CXX) $(OBJ_CXXFLAGS) $(OSX_LDFLAGS_i686) -o milkyway_$(APP_VERSION)_i686-apple-darwin $(APP_OBJS) $(SEARCH_OBJS) $(UTIL_OBJS) -lboinc -lboinc_api

osx_i686_gpu: OBJ_CXX = $(OSX_CXX_i686)
osx_i686_gpu: OBJ_CXXFLAGS = $(GPU_FLAGS) $(OSX_CXXFLAGS_i686)
osx_i686_gpu: $(GPU_APP_OBJS) $(APP_OBJS) $(SEARCH_OBJS) $(UTIL_OBJS) $(OSX_LIBBOINC_DIR)/libboinc_api.a $(OSX_LIBBOINC_DIR)/libboinc.a
	$(OBJ_CXX) $(OBJ_CXXFLAGS) $(LDFLAGS_x86_64) -o milkywayGPU_$(APP_VERSION)_i686-apple-darwin $(APP_OBJS) $(SEARCH_OBJS) $(UTIL_OBJS) $(GPU_APP_OBJS) -lboinc -lboinc_api -lcudart -lcutil

osx_x86_64: OBJ_CXX = $(OSX_CXX_x86_64)
osx_x86_64: OBJ_CXXFLAGS = $(CPU_FLAGS) $(OSX_CXXFLAGS_x86_64)
osx_x86_64: $(APP_OBJS) $(SEARCH_OBJS) $(UTIL_OBJS) $(OSX_LIBBOINC_DIR)/libboinc_api.a $(OSX_LIBBOINC_DIR)/libboinc.a
	$(OBJ_CXX) $(OBJ_CXXFLAGS) $(OSX_LDFLAGS_x86_64) -o milkyway_$(APP_VERSION)_x86_64-apple-darwin $(APP_OBJS) $(SEARCH_OBJS) $(UTIL_OBJS) -lboinc -lboinc_api 


linux_i686: OBJ_CXX = $(LINUX_CXX_i686)
linux_i686: OBJ_CXXFLAGS = $(CPU_FLAGS) $(LINUX_CXXFLAGS_i686)
linux_i686: $(APP_OBJS) $(SEARCH_OBJS) $(UTIL_OBJS) $(LINUX_LIBBOINC_DIR)/libboinc_api.a $(LINUX_LIBBOINC_DIR)/libboinc.a
	$(OBJ_CXX) $(OBJ_CXXFLAGS) $(LINUX_LDFLAGS_i686) -Wl --export_dynamic -o milkyway_$(APP_VERSION)_i686-pc-linux-gnu $(APP_OBJS) $(SEARCH_OBJS) $(UTIL_OBJS) -lm -lboinc_api -lboinc -pthread

linux_i686_gpu: OBJ_CXX = $(LINUX_CXX_i686)
linux_i686_gpu: OBJ_CXXFLAGS = $(GPU_FLAGS) $(LINUX_CXXFLAGS_i686)
linux_i686_gpu: $(GPU_APP_OBJS) $(APP_OBJS) $(SEARCH_OBJS) $(UTIL_OBJS) $(LINUX_LIBBOINC_DIR)/libboinc_api.a $(LINUX_LIBBOINC_DIR)/libboinc.a
	$(OBJ_CXX) $(OBJ_CXXFLAGS) $(LDFLAGS_x86_64) -o milkywayGPU_$(APP_VERSION)_i686-pc-linux-gnu $(APP_OBJS) $(SEARCH_OBJS) $(UTIL_OBJS) $(GPU_APP_OBJS) -lboinc -lboinc_api -lcudart -lcutil

linux_x86_64: OBJ_CXX = $(LINUX_CXX_x86_64)
linux_x86_64: OBJ_CXXFLAGS = $(CPU_FLAGS) $(LINUX_CXXFLAGS_x86_64)
linux_x86_64: $(APP_OBJS) $(SEARCH_OBJS) $(UTIL_OBJS) $(LINUX_LIBBOINC_DIR)/libboinc_api.a $(LINUX_LIBBOINC_DIR)/libboinc.a
	$(OBJ_CXX) $(OBJ_CXXFLAGS) $(LINUX_LDFLAGS_x86_64) -Wl -static -static-libgcc -o milkyway_$(APP_VERSION)_x86_64-pc-linux-gnu $(APP_OBJS) $(SEARCH_OBJS) $(UTIL_OBJS) -lm -lboinc_api -lboinc -pthread

linux_x86_64_gpu: OBJ_CXX = $(LINUX_CXX_x86_64)
linux_x86_64_gpu: OBJ_CXXFLAGS = $(GPU_FLAGS) $(LINUX_CXXFLAGS_x86_64)
linux_x86_64_gpu: $(GPU_APP_OBJS) $(APP_OBJS) $(SEARCH_OBJS) $(UTIL_OBJS) 
	$(OBJ_CXX) $(OBJ_CXXFLAGS) $(LDFLAGS_x86_64) -o milkywayGPU_$(APP_VERSION)_x86_64-pc-linux-gnu $(APP_OBJS) $(SEARCH_OBJS) $(UTIL_OBJS) $(GPU_APP_OBJS) -lboinc_api -lboinc -lcudart -lcutil


.C.o:
	$(OBJ_CXX) $(OBJ_CXXFLAGS) $(INC) -Wall -c $< -o $@
.c.o:
	$(OBJ_CXX) $(OBJ_CXXFLAGS) $(INC) -Wall -c $< -o $@


clean:
	rm -f $(SEARCH_OBJS) $(UTIL_OBJS) $(APP_OBJS) $(GPU_APP_OBJS) evaluation_gpu.linkinfo evaluation_gpu2.linkinfo evaluation_gpu3.linkinfo

clean_all:
	rm -f $(PROGS) $(SEARCH_OBJS) $(UTIL_OBJS) $(APP_OBJS)

################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../separation/src/separation_graphics.cc 

CPP_SRCS += \
../separation/src/cl_compile_flags.cpp 

C_SRCS += \
../separation/src/calculated_constants.c \
../separation/src/coordinates.c \
../separation/src/evaluation.c \
../separation/src/evaluation_state.c \
../separation/src/gauss_legendre.c \
../separation/src/il_kernel_1.c \
../separation/src/il_kernel_2.c \
../separation/src/il_kernel_3.c \
../separation/src/il_kernel_4.c \
../separation/src/integrals.c \
../separation/src/io_util.c \
../separation/src/likelihood.c \
../separation/src/parameters.c \
../separation/src/probabilities.c \
../separation/src/probabilities_avx.c \
../separation/src/probabilities_dispatch.c \
../separation/src/probabilities_sse2.c \
../separation/src/r_points.c \
../separation/src/replace_amd_il.c \
../separation/src/run_cl.c \
../separation/src/separation_binaries.c \
../separation/src/separation_cl_buffers.c \
../separation/src/separation_lua.c \
../separation/src/separation_main.c \
../separation/src/separation_utils.c \
../separation/src/setup_cl.c \
../separation/src/star_points.c 

CC_DEPS += \
./separation/src/separation_graphics.d 

OBJS += \
./separation/src/calculated_constants.o \
./separation/src/cl_compile_flags.o \
./separation/src/coordinates.o \
./separation/src/evaluation.o \
./separation/src/evaluation_state.o \
./separation/src/gauss_legendre.o \
./separation/src/il_kernel_1.o \
./separation/src/il_kernel_2.o \
./separation/src/il_kernel_3.o \
./separation/src/il_kernel_4.o \
./separation/src/integrals.o \
./separation/src/io_util.o \
./separation/src/likelihood.o \
./separation/src/parameters.o \
./separation/src/probabilities.o \
./separation/src/probabilities_avx.o \
./separation/src/probabilities_dispatch.o \
./separation/src/probabilities_sse2.o \
./separation/src/r_points.o \
./separation/src/replace_amd_il.o \
./separation/src/run_cl.o \
./separation/src/separation_binaries.o \
./separation/src/separation_cl_buffers.o \
./separation/src/separation_graphics.o \
./separation/src/separation_lua.o \
./separation/src/separation_main.o \
./separation/src/separation_utils.o \
./separation/src/setup_cl.o \
./separation/src/star_points.o 

CPP_DEPS += \
./separation/src/cl_compile_flags.d 

C_DEPS += \
./separation/src/calculated_constants.d \
./separation/src/coordinates.d \
./separation/src/evaluation.d \
./separation/src/evaluation_state.d \
./separation/src/gauss_legendre.d \
./separation/src/il_kernel_1.d \
./separation/src/il_kernel_2.d \
./separation/src/il_kernel_3.d \
./separation/src/il_kernel_4.d \
./separation/src/integrals.d \
./separation/src/io_util.d \
./separation/src/likelihood.d \
./separation/src/parameters.d \
./separation/src/probabilities.d \
./separation/src/probabilities_avx.d \
./separation/src/probabilities_dispatch.d \
./separation/src/probabilities_sse2.d \
./separation/src/r_points.d \
./separation/src/replace_amd_il.d \
./separation/src/run_cl.d \
./separation/src/separation_binaries.d \
./separation/src/separation_cl_buffers.d \
./separation/src/separation_lua.d \
./separation/src/separation_main.d \
./separation/src/separation_utils.d \
./separation/src/setup_cl.d \
./separation/src/star_points.d 


# Each subdirectory must supply rules for building sources it contributes
separation/src/%.o: ../separation/src/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

separation/src/%.o: ../separation/src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

separation/src/%.o: ../separation/src/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



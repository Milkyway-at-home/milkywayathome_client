################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../glm/test/core/core_func_common.cpp \
../glm/test/core/core_func_exponential.cpp \
../glm/test/core/core_func_geometric.cpp \
../glm/test/core/core_func_integer.cpp \
../glm/test/core/core_func_matrix.cpp \
../glm/test/core/core_func_noise.cpp \
../glm/test/core/core_func_packing.cpp \
../glm/test/core/core_func_swizzle.cpp \
../glm/test/core/core_func_trigonometric.cpp \
../glm/test/core/core_func_vector_relational.cpp \
../glm/test/core/core_setup_message.cpp \
../glm/test/core/core_setup_precision.cpp \
../glm/test/core/core_type_float.cpp \
../glm/test/core/core_type_half.cpp \
../glm/test/core/core_type_int.cpp \
../glm/test/core/core_type_length.cpp \
../glm/test/core/core_type_mat2x2.cpp \
../glm/test/core/core_type_mat2x3.cpp \
../glm/test/core/core_type_mat2x4.cpp \
../glm/test/core/core_type_mat3x2.cpp \
../glm/test/core/core_type_mat3x3.cpp \
../glm/test/core/core_type_mat3x4.cpp \
../glm/test/core/core_type_mat4x2.cpp \
../glm/test/core/core_type_mat4x3.cpp \
../glm/test/core/core_type_mat4x4.cpp \
../glm/test/core/core_type_vec1.cpp \
../glm/test/core/core_type_vec2.cpp \
../glm/test/core/core_type_vec3.cpp \
../glm/test/core/core_type_vec4.cpp 

OBJS += \
./glm/test/core/core_func_common.o \
./glm/test/core/core_func_exponential.o \
./glm/test/core/core_func_geometric.o \
./glm/test/core/core_func_integer.o \
./glm/test/core/core_func_matrix.o \
./glm/test/core/core_func_noise.o \
./glm/test/core/core_func_packing.o \
./glm/test/core/core_func_swizzle.o \
./glm/test/core/core_func_trigonometric.o \
./glm/test/core/core_func_vector_relational.o \
./glm/test/core/core_setup_message.o \
./glm/test/core/core_setup_precision.o \
./glm/test/core/core_type_float.o \
./glm/test/core/core_type_half.o \
./glm/test/core/core_type_int.o \
./glm/test/core/core_type_length.o \
./glm/test/core/core_type_mat2x2.o \
./glm/test/core/core_type_mat2x3.o \
./glm/test/core/core_type_mat2x4.o \
./glm/test/core/core_type_mat3x2.o \
./glm/test/core/core_type_mat3x3.o \
./glm/test/core/core_type_mat3x4.o \
./glm/test/core/core_type_mat4x2.o \
./glm/test/core/core_type_mat4x3.o \
./glm/test/core/core_type_mat4x4.o \
./glm/test/core/core_type_vec1.o \
./glm/test/core/core_type_vec2.o \
./glm/test/core/core_type_vec3.o \
./glm/test/core/core_type_vec4.o 

CPP_DEPS += \
./glm/test/core/core_func_common.d \
./glm/test/core/core_func_exponential.d \
./glm/test/core/core_func_geometric.d \
./glm/test/core/core_func_integer.d \
./glm/test/core/core_func_matrix.d \
./glm/test/core/core_func_noise.d \
./glm/test/core/core_func_packing.d \
./glm/test/core/core_func_swizzle.d \
./glm/test/core/core_func_trigonometric.d \
./glm/test/core/core_func_vector_relational.d \
./glm/test/core/core_setup_message.d \
./glm/test/core/core_setup_precision.d \
./glm/test/core/core_type_float.d \
./glm/test/core/core_type_half.d \
./glm/test/core/core_type_int.d \
./glm/test/core/core_type_length.d \
./glm/test/core/core_type_mat2x2.d \
./glm/test/core/core_type_mat2x3.d \
./glm/test/core/core_type_mat2x4.d \
./glm/test/core/core_type_mat3x2.d \
./glm/test/core/core_type_mat3x3.d \
./glm/test/core/core_type_mat3x4.d \
./glm/test/core/core_type_mat4x2.d \
./glm/test/core/core_type_mat4x3.d \
./glm/test/core/core_type_mat4x4.d \
./glm/test/core/core_type_vec1.d \
./glm/test/core/core_type_vec2.d \
./glm/test/core/core_type_vec3.d \
./glm/test/core/core_type_vec4.d 


# Each subdirectory must supply rules for building sources it contributes
glm/test/core/%.o: ../glm/test/core/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../glm/test/gtx/gtx_bit.cpp \
../glm/test/gtx/gtx_gradient_paint.cpp \
../glm/test/gtx/gtx_integer.cpp \
../glm/test/gtx/gtx_matrix_query.cpp \
../glm/test/gtx/gtx_noise.cpp \
../glm/test/gtx/gtx_quaternion.cpp \
../glm/test/gtx/gtx_random.cpp \
../glm/test/gtx/gtx_rotate_vector.cpp \
../glm/test/gtx/gtx_simd_mat4.cpp \
../glm/test/gtx/gtx_simd_vec4.cpp \
../glm/test/gtx/gtx_string_cast.cpp \
../glm/test/gtx/gtx_ulp.cpp \
../glm/test/gtx/gtx_vector_angle.cpp \
../glm/test/gtx/gtx_vector_query.cpp 

OBJS += \
./glm/test/gtx/gtx_bit.o \
./glm/test/gtx/gtx_gradient_paint.o \
./glm/test/gtx/gtx_integer.o \
./glm/test/gtx/gtx_matrix_query.o \
./glm/test/gtx/gtx_noise.o \
./glm/test/gtx/gtx_quaternion.o \
./glm/test/gtx/gtx_random.o \
./glm/test/gtx/gtx_rotate_vector.o \
./glm/test/gtx/gtx_simd_mat4.o \
./glm/test/gtx/gtx_simd_vec4.o \
./glm/test/gtx/gtx_string_cast.o \
./glm/test/gtx/gtx_ulp.o \
./glm/test/gtx/gtx_vector_angle.o \
./glm/test/gtx/gtx_vector_query.o 

CPP_DEPS += \
./glm/test/gtx/gtx_bit.d \
./glm/test/gtx/gtx_gradient_paint.d \
./glm/test/gtx/gtx_integer.d \
./glm/test/gtx/gtx_matrix_query.d \
./glm/test/gtx/gtx_noise.d \
./glm/test/gtx/gtx_quaternion.d \
./glm/test/gtx/gtx_random.d \
./glm/test/gtx/gtx_rotate_vector.d \
./glm/test/gtx/gtx_simd_mat4.d \
./glm/test/gtx/gtx_simd_vec4.d \
./glm/test/gtx/gtx_string_cast.d \
./glm/test/gtx/gtx_ulp.d \
./glm/test/gtx/gtx_vector_angle.d \
./glm/test/gtx/gtx_vector_query.d 


# Each subdirectory must supply rules for building sources it contributes
glm/test/gtx/%.o: ../glm/test/gtx/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../glm/test/gtc/gtc_half_float.cpp \
../glm/test/gtc/gtc_matrix_access.cpp \
../glm/test/gtc/gtc_matrix_integer.cpp \
../glm/test/gtc/gtc_matrix_inverse.cpp \
../glm/test/gtc/gtc_matrix_transform.cpp \
../glm/test/gtc/gtc_noise.cpp \
../glm/test/gtc/gtc_quaternion.cpp \
../glm/test/gtc/gtc_random.cpp \
../glm/test/gtc/gtc_swizzle.cpp \
../glm/test/gtc/gtc_type_precision.cpp \
../glm/test/gtc/gtc_type_ptr.cpp 

OBJS += \
./glm/test/gtc/gtc_half_float.o \
./glm/test/gtc/gtc_matrix_access.o \
./glm/test/gtc/gtc_matrix_integer.o \
./glm/test/gtc/gtc_matrix_inverse.o \
./glm/test/gtc/gtc_matrix_transform.o \
./glm/test/gtc/gtc_noise.o \
./glm/test/gtc/gtc_quaternion.o \
./glm/test/gtc/gtc_random.o \
./glm/test/gtc/gtc_swizzle.o \
./glm/test/gtc/gtc_type_precision.o \
./glm/test/gtc/gtc_type_ptr.o 

CPP_DEPS += \
./glm/test/gtc/gtc_half_float.d \
./glm/test/gtc/gtc_matrix_access.d \
./glm/test/gtc/gtc_matrix_integer.d \
./glm/test/gtc/gtc_matrix_inverse.d \
./glm/test/gtc/gtc_matrix_transform.d \
./glm/test/gtc/gtc_noise.d \
./glm/test/gtc/gtc_quaternion.d \
./glm/test/gtc/gtc_random.d \
./glm/test/gtc/gtc_swizzle.d \
./glm/test/gtc/gtc_type_precision.d \
./glm/test/gtc/gtc_type_ptr.d 


# Each subdirectory must supply rules for building sources it contributes
glm/test/gtc/%.o: ../glm/test/gtc/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../glm/doc/src/dummy.cpp 

OBJS += \
./glm/doc/src/dummy.o 

CPP_DEPS += \
./glm/doc/src/dummy.d 


# Each subdirectory must supply rules for building sources it contributes
glm/doc/src/%.o: ../glm/doc/src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../glm/util/glm_core.cpp 

OBJS += \
./glm/util/glm_core.o 

CPP_DEPS += \
./glm/util/glm_core.d 


# Each subdirectory must supply rules for building sources it contributes
glm/util/%.o: ../glm/util/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



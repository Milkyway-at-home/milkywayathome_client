################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../glm/test/external/gli/core/dummy.cpp 

OBJS += \
./glm/test/external/gli/core/dummy.o 

CPP_DEPS += \
./glm/test/external/gli/core/dummy.d 


# Each subdirectory must supply rules for building sources it contributes
glm/test/external/gli/core/%.o: ../glm/test/external/gli/core/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



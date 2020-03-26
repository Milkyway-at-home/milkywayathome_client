################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../tools/separation_il_kernelgen/il_kernelgen.cpp 

OBJS += \
./tools/separation_il_kernelgen/il_kernelgen.o 

CPP_DEPS += \
./tools/separation_il_kernelgen/il_kernelgen.d 


# Each subdirectory must supply rules for building sources it contributes
tools/separation_il_kernelgen/%.o: ../tools/separation_il_kernelgen/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



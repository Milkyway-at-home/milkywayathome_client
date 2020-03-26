################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CXX_SRCS += \
../mwtao/mle.cxx 

CXX_DEPS += \
./mwtao/mle.d 

OBJS += \
./mwtao/mle.o 


# Each subdirectory must supply rules for building sources it contributes
mwtao/%.o: ../mwtao/%.cxx
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



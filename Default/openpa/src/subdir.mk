################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../openpa/src/opa_primitives.c \
../openpa/src/opa_queue.c 

OBJS += \
./openpa/src/opa_primitives.o \
./openpa/src/opa_queue.o 

C_DEPS += \
./openpa/src/opa_primitives.d \
./openpa/src/opa_queue.d 


# Each subdirectory must supply rules for building sources it contributes
openpa/src/%.o: ../openpa/src/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../openpa/test/sanity.c \
../openpa/test/test_barriers.c \
../openpa/test/test_primitives.c \
../openpa/test/test_queue.c 

OBJS += \
./openpa/test/sanity.o \
./openpa/test/test_barriers.o \
./openpa/test/test_primitives.o \
./openpa/test/test_queue.o 

C_DEPS += \
./openpa/test/sanity.d \
./openpa/test/test_barriers.d \
./openpa/test/test_primitives.d \
./openpa/test/test_queue.d 


# Each subdirectory must supply rules for building sources it contributes
openpa/test/%.o: ../openpa/test/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



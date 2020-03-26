################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../popt/src/lookup3.c \
../popt/src/popt.c \
../popt/src/poptconfig.c \
../popt/src/popthelp.c \
../popt/src/poptint.c \
../popt/src/poptparse.c 

OBJS += \
./popt/src/lookup3.o \
./popt/src/popt.o \
./popt/src/poptconfig.o \
./popt/src/popthelp.o \
./popt/src/poptint.o \
./popt/src/poptparse.o 

C_DEPS += \
./popt/src/lookup3.d \
./popt/src/popt.d \
./popt/src/poptconfig.d \
./popt/src/popthelp.d \
./popt/src/poptint.d \
./popt/src/poptparse.d 


# Each subdirectory must supply rules for building sources it contributes
popt/src/%.o: ../popt/src/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



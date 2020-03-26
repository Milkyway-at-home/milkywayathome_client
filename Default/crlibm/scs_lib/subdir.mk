################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../crlibm/scs_lib/addition_scs.c \
../crlibm/scs_lib/division_scs.c \
../crlibm/scs_lib/double2scs.c \
../crlibm/scs_lib/multiplication_scs.c \
../crlibm/scs_lib/poly_fct.c \
../crlibm/scs_lib/print_scs.c \
../crlibm/scs_lib/rand_scs.c \
../crlibm/scs_lib/scs2double.c \
../crlibm/scs_lib/scs2mpf.c \
../crlibm/scs_lib/scs2mpfr.c \
../crlibm/scs_lib/scs_private.c \
../crlibm/scs_lib/zero_scs.c 

OBJS += \
./crlibm/scs_lib/addition_scs.o \
./crlibm/scs_lib/division_scs.o \
./crlibm/scs_lib/double2scs.o \
./crlibm/scs_lib/multiplication_scs.o \
./crlibm/scs_lib/poly_fct.o \
./crlibm/scs_lib/print_scs.o \
./crlibm/scs_lib/rand_scs.o \
./crlibm/scs_lib/scs2double.o \
./crlibm/scs_lib/scs2mpf.o \
./crlibm/scs_lib/scs2mpfr.o \
./crlibm/scs_lib/scs_private.o \
./crlibm/scs_lib/zero_scs.o 

C_DEPS += \
./crlibm/scs_lib/addition_scs.d \
./crlibm/scs_lib/division_scs.d \
./crlibm/scs_lib/double2scs.d \
./crlibm/scs_lib/multiplication_scs.d \
./crlibm/scs_lib/poly_fct.d \
./crlibm/scs_lib/print_scs.d \
./crlibm/scs_lib/rand_scs.d \
./crlibm/scs_lib/scs2double.d \
./crlibm/scs_lib/scs2mpf.d \
./crlibm/scs_lib/scs2mpfr.d \
./crlibm/scs_lib/scs_private.d \
./crlibm/scs_lib/zero_scs.d 


# Each subdirectory must supply rules for building sources it contributes
crlibm/scs_lib/%.o: ../crlibm/scs_lib/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



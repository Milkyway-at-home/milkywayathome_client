################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../crlibm/acos-td.c \
../crlibm/asin-td.c \
../crlibm/asincos.c \
../crlibm/atan_accurate.c \
../crlibm/atan_fast.c \
../crlibm/crlibm_private.c \
../crlibm/csh_fast.c \
../crlibm/exp-td-standalone.c \
../crlibm/exp-td.c \
../crlibm/expm1-standalone.c \
../crlibm/expm1.c \
../crlibm/log-de.c \
../crlibm/log-de2.c \
../crlibm/log-td.c \
../crlibm/log.c \
../crlibm/log10-td.c \
../crlibm/log1p.c \
../crlibm/log2-td.c \
../crlibm/log2_accurate.c \
../crlibm/log_accurate.c \
../crlibm/log_fast.c \
../crlibm/pow.c \
../crlibm/rem_pio2_accurate.c \
../crlibm/trigo_accurate.c \
../crlibm/trigo_fast.c \
../crlibm/trigpi.c \
../crlibm/triple-double.c 

OBJS += \
./crlibm/acos-td.o \
./crlibm/asin-td.o \
./crlibm/asincos.o \
./crlibm/atan_accurate.o \
./crlibm/atan_fast.o \
./crlibm/crlibm_private.o \
./crlibm/csh_fast.o \
./crlibm/exp-td-standalone.o \
./crlibm/exp-td.o \
./crlibm/expm1-standalone.o \
./crlibm/expm1.o \
./crlibm/log-de.o \
./crlibm/log-de2.o \
./crlibm/log-td.o \
./crlibm/log.o \
./crlibm/log10-td.o \
./crlibm/log1p.o \
./crlibm/log2-td.o \
./crlibm/log2_accurate.o \
./crlibm/log_accurate.o \
./crlibm/log_fast.o \
./crlibm/pow.o \
./crlibm/rem_pio2_accurate.o \
./crlibm/trigo_accurate.o \
./crlibm/trigo_fast.o \
./crlibm/trigpi.o \
./crlibm/triple-double.o 

C_DEPS += \
./crlibm/acos-td.d \
./crlibm/asin-td.d \
./crlibm/asincos.d \
./crlibm/atan_accurate.d \
./crlibm/atan_fast.d \
./crlibm/crlibm_private.d \
./crlibm/csh_fast.d \
./crlibm/exp-td-standalone.d \
./crlibm/exp-td.d \
./crlibm/expm1-standalone.d \
./crlibm/expm1.d \
./crlibm/log-de.d \
./crlibm/log-de2.d \
./crlibm/log-td.d \
./crlibm/log.d \
./crlibm/log10-td.d \
./crlibm/log1p.d \
./crlibm/log2-td.d \
./crlibm/log2_accurate.d \
./crlibm/log_accurate.d \
./crlibm/log_fast.d \
./crlibm/pow.d \
./crlibm/rem_pio2_accurate.d \
./crlibm/trigo_accurate.d \
./crlibm/trigo_fast.d \
./crlibm/trigpi.d \
./crlibm/triple-double.d 


# Each subdirectory must supply rules for building sources it contributes
crlibm/%.o: ../crlibm/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../crlibm/tests/blind_test.c \
../crlibm/tests/filterlists.c \
../crlibm/tests/generate_test_vectors.c \
../crlibm/tests/interval_explog.c \
../crlibm/tests/soak_test-interval.c \
../crlibm/tests/soak_test.c \
../crlibm/tests/soak_test_trigpi.c \
../crlibm/tests/tanbug.c \
../crlibm/tests/test_common.c \
../crlibm/tests/test_perf-interval.c \
../crlibm/tests/test_perf.c \
../crlibm/tests/test_val.c 

OBJS += \
./crlibm/tests/blind_test.o \
./crlibm/tests/filterlists.o \
./crlibm/tests/generate_test_vectors.o \
./crlibm/tests/interval_explog.o \
./crlibm/tests/soak_test-interval.o \
./crlibm/tests/soak_test.o \
./crlibm/tests/soak_test_trigpi.o \
./crlibm/tests/tanbug.o \
./crlibm/tests/test_common.o \
./crlibm/tests/test_perf-interval.o \
./crlibm/tests/test_perf.o \
./crlibm/tests/test_val.o 

C_DEPS += \
./crlibm/tests/blind_test.d \
./crlibm/tests/filterlists.d \
./crlibm/tests/generate_test_vectors.d \
./crlibm/tests/interval_explog.d \
./crlibm/tests/soak_test-interval.d \
./crlibm/tests/soak_test.d \
./crlibm/tests/soak_test_trigpi.d \
./crlibm/tests/tanbug.d \
./crlibm/tests/test_common.d \
./crlibm/tests/test_perf-interval.d \
./crlibm/tests/test_perf.d \
./crlibm/tests/test_val.d 


# Each subdirectory must supply rules for building sources it contributes
crlibm/tests/%.o: ../crlibm/tests/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



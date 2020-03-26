################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../boinc/boinc/zip/zip/__p___mb_cur_max.c \
../boinc/boinc/zip/zip/deflate.c \
../boinc/boinc/zip/zip/trees.c \
../boinc/boinc/zip/zip/util.c \
../boinc/boinc/zip/zip/z_fileio.c \
../boinc/boinc/zip/zip/z_globals.c \
../boinc/boinc/zip/zip/zip.c \
../boinc/boinc/zip/zip/zipfile.c \
../boinc/boinc/zip/zip/zipup.c 

OBJS += \
./boinc/boinc/zip/zip/__p___mb_cur_max.o \
./boinc/boinc/zip/zip/deflate.o \
./boinc/boinc/zip/zip/trees.o \
./boinc/boinc/zip/zip/util.o \
./boinc/boinc/zip/zip/z_fileio.o \
./boinc/boinc/zip/zip/z_globals.o \
./boinc/boinc/zip/zip/zip.o \
./boinc/boinc/zip/zip/zipfile.o \
./boinc/boinc/zip/zip/zipup.o 

C_DEPS += \
./boinc/boinc/zip/zip/__p___mb_cur_max.d \
./boinc/boinc/zip/zip/deflate.d \
./boinc/boinc/zip/zip/trees.d \
./boinc/boinc/zip/zip/util.d \
./boinc/boinc/zip/zip/z_fileio.d \
./boinc/boinc/zip/zip/z_globals.d \
./boinc/boinc/zip/zip/zip.d \
./boinc/boinc/zip/zip/zipfile.d \
./boinc/boinc/zip/zip/zipup.d 


# Each subdirectory must supply rules for building sources it contributes
boinc/boinc/zip/zip/%.o: ../boinc/boinc/zip/zip/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



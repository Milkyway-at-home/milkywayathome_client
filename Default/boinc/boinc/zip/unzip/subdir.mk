################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../boinc/boinc/zip/unzip/api.c \
../boinc/boinc/zip/unzip/apihelp.c \
../boinc/boinc/zip/unzip/crc32.c \
../boinc/boinc/zip/unzip/crypt.c \
../boinc/boinc/zip/unzip/explode.c \
../boinc/boinc/zip/unzip/extract.c \
../boinc/boinc/zip/unzip/fileio.c \
../boinc/boinc/zip/unzip/globals.c \
../boinc/boinc/zip/unzip/inflate.c \
../boinc/boinc/zip/unzip/list.c \
../boinc/boinc/zip/unzip/match.c \
../boinc/boinc/zip/unzip/process.c \
../boinc/boinc/zip/unzip/ttyio.c \
../boinc/boinc/zip/unzip/unreduce.c \
../boinc/boinc/zip/unzip/unshrink.c \
../boinc/boinc/zip/unzip/unzip.c \
../boinc/boinc/zip/unzip/zipinfo.c 

OBJS += \
./boinc/boinc/zip/unzip/api.o \
./boinc/boinc/zip/unzip/apihelp.o \
./boinc/boinc/zip/unzip/crc32.o \
./boinc/boinc/zip/unzip/crypt.o \
./boinc/boinc/zip/unzip/explode.o \
./boinc/boinc/zip/unzip/extract.o \
./boinc/boinc/zip/unzip/fileio.o \
./boinc/boinc/zip/unzip/globals.o \
./boinc/boinc/zip/unzip/inflate.o \
./boinc/boinc/zip/unzip/list.o \
./boinc/boinc/zip/unzip/match.o \
./boinc/boinc/zip/unzip/process.o \
./boinc/boinc/zip/unzip/ttyio.o \
./boinc/boinc/zip/unzip/unreduce.o \
./boinc/boinc/zip/unzip/unshrink.o \
./boinc/boinc/zip/unzip/unzip.o \
./boinc/boinc/zip/unzip/zipinfo.o 

C_DEPS += \
./boinc/boinc/zip/unzip/api.d \
./boinc/boinc/zip/unzip/apihelp.d \
./boinc/boinc/zip/unzip/crc32.d \
./boinc/boinc/zip/unzip/crypt.d \
./boinc/boinc/zip/unzip/explode.d \
./boinc/boinc/zip/unzip/extract.d \
./boinc/boinc/zip/unzip/fileio.d \
./boinc/boinc/zip/unzip/globals.d \
./boinc/boinc/zip/unzip/inflate.d \
./boinc/boinc/zip/unzip/list.d \
./boinc/boinc/zip/unzip/match.d \
./boinc/boinc/zip/unzip/process.d \
./boinc/boinc/zip/unzip/ttyio.d \
./boinc/boinc/zip/unzip/unreduce.d \
./boinc/boinc/zip/unzip/unshrink.d \
./boinc/boinc/zip/unzip/unzip.d \
./boinc/boinc/zip/unzip/zipinfo.d 


# Each subdirectory must supply rules for building sources it contributes
boinc/boinc/zip/unzip/%.o: ../boinc/boinc/zip/unzip/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



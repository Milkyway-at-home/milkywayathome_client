################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../boinc/boinc/samples/glut/glut_roman.c \
../boinc/boinc/samples/glut/glut_stroke.c \
../boinc/boinc/samples/glut/glut_swidth.c \
../boinc/boinc/samples/glut/win32_glx.c \
../boinc/boinc/samples/glut/win32_util.c \
../boinc/boinc/samples/glut/win32_x11.c 

OBJS += \
./boinc/boinc/samples/glut/glut_roman.o \
./boinc/boinc/samples/glut/glut_stroke.o \
./boinc/boinc/samples/glut/glut_swidth.o \
./boinc/boinc/samples/glut/win32_glx.o \
./boinc/boinc/samples/glut/win32_util.o \
./boinc/boinc/samples/glut/win32_x11.o 

C_DEPS += \
./boinc/boinc/samples/glut/glut_roman.d \
./boinc/boinc/samples/glut/glut_stroke.d \
./boinc/boinc/samples/glut/glut_swidth.d \
./boinc/boinc/samples/glut/win32_glx.d \
./boinc/boinc/samples/glut/win32_util.d \
./boinc/boinc/samples/glut/win32_x11.d 


# Each subdirectory must supply rules for building sources it contributes
boinc/boinc/samples/glut/%.o: ../boinc/boinc/samples/glut/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../glfw/tests/accuracy.c \
../glfw/tests/clipboard.c \
../glfw/tests/defaults.c \
../glfw/tests/events.c \
../glfw/tests/fsaa.c \
../glfw/tests/gamma.c \
../glfw/tests/glfwinfo.c \
../glfw/tests/iconify.c \
../glfw/tests/joysticks.c \
../glfw/tests/modes.c \
../glfw/tests/peter.c \
../glfw/tests/reopen.c \
../glfw/tests/sharing.c \
../glfw/tests/tearing.c \
../glfw/tests/threads.c \
../glfw/tests/title.c \
../glfw/tests/windows.c 

OBJS += \
./glfw/tests/accuracy.o \
./glfw/tests/clipboard.o \
./glfw/tests/defaults.o \
./glfw/tests/events.o \
./glfw/tests/fsaa.o \
./glfw/tests/gamma.o \
./glfw/tests/glfwinfo.o \
./glfw/tests/iconify.o \
./glfw/tests/joysticks.o \
./glfw/tests/modes.o \
./glfw/tests/peter.o \
./glfw/tests/reopen.o \
./glfw/tests/sharing.o \
./glfw/tests/tearing.o \
./glfw/tests/threads.o \
./glfw/tests/title.o \
./glfw/tests/windows.o 

C_DEPS += \
./glfw/tests/accuracy.d \
./glfw/tests/clipboard.d \
./glfw/tests/defaults.d \
./glfw/tests/events.d \
./glfw/tests/fsaa.d \
./glfw/tests/gamma.d \
./glfw/tests/glfwinfo.d \
./glfw/tests/iconify.d \
./glfw/tests/joysticks.d \
./glfw/tests/modes.d \
./glfw/tests/peter.d \
./glfw/tests/reopen.d \
./glfw/tests/sharing.d \
./glfw/tests/tearing.d \
./glfw/tests/threads.d \
./glfw/tests/title.d \
./glfw/tests/windows.d 


# Each subdirectory must supply rules for building sources it contributes
glfw/tests/%.o: ../glfw/tests/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



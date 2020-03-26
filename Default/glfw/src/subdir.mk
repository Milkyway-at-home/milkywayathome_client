################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../glfw/src/clipboard.c \
../glfw/src/cocoa_gamma.c \
../glfw/src/cocoa_time.c \
../glfw/src/context.c \
../glfw/src/egl_context.c \
../glfw/src/gamma.c \
../glfw/src/glx_context.c \
../glfw/src/init.c \
../glfw/src/input.c \
../glfw/src/joystick.c \
../glfw/src/monitor.c \
../glfw/src/time.c \
../glfw/src/wgl_context.c \
../glfw/src/win32_clipboard.c \
../glfw/src/win32_gamma.c \
../glfw/src/win32_init.c \
../glfw/src/win32_joystick.c \
../glfw/src/win32_monitor.c \
../glfw/src/win32_time.c \
../glfw/src/win32_window.c \
../glfw/src/window.c \
../glfw/src/x11_clipboard.c \
../glfw/src/x11_gamma.c \
../glfw/src/x11_init.c \
../glfw/src/x11_joystick.c \
../glfw/src/x11_monitor.c \
../glfw/src/x11_time.c \
../glfw/src/x11_unicode.c \
../glfw/src/x11_window.c 

OBJS += \
./glfw/src/clipboard.o \
./glfw/src/cocoa_gamma.o \
./glfw/src/cocoa_time.o \
./glfw/src/context.o \
./glfw/src/egl_context.o \
./glfw/src/gamma.o \
./glfw/src/glx_context.o \
./glfw/src/init.o \
./glfw/src/input.o \
./glfw/src/joystick.o \
./glfw/src/monitor.o \
./glfw/src/time.o \
./glfw/src/wgl_context.o \
./glfw/src/win32_clipboard.o \
./glfw/src/win32_gamma.o \
./glfw/src/win32_init.o \
./glfw/src/win32_joystick.o \
./glfw/src/win32_monitor.o \
./glfw/src/win32_time.o \
./glfw/src/win32_window.o \
./glfw/src/window.o \
./glfw/src/x11_clipboard.o \
./glfw/src/x11_gamma.o \
./glfw/src/x11_init.o \
./glfw/src/x11_joystick.o \
./glfw/src/x11_monitor.o \
./glfw/src/x11_time.o \
./glfw/src/x11_unicode.o \
./glfw/src/x11_window.o 

C_DEPS += \
./glfw/src/clipboard.d \
./glfw/src/cocoa_gamma.d \
./glfw/src/cocoa_time.d \
./glfw/src/context.d \
./glfw/src/egl_context.d \
./glfw/src/gamma.d \
./glfw/src/glx_context.d \
./glfw/src/init.d \
./glfw/src/input.d \
./glfw/src/joystick.d \
./glfw/src/monitor.d \
./glfw/src/time.d \
./glfw/src/wgl_context.d \
./glfw/src/win32_clipboard.d \
./glfw/src/win32_gamma.d \
./glfw/src/win32_init.d \
./glfw/src/win32_joystick.d \
./glfw/src/win32_monitor.d \
./glfw/src/win32_time.d \
./glfw/src/win32_window.d \
./glfw/src/window.d \
./glfw/src/x11_clipboard.d \
./glfw/src/x11_gamma.d \
./glfw/src/x11_init.d \
./glfw/src/x11_joystick.d \
./glfw/src/x11_monitor.d \
./glfw/src/x11_time.d \
./glfw/src/x11_unicode.d \
./glfw/src/x11_window.d 


# Each subdirectory must supply rules for building sources it contributes
glfw/src/%.o: ../glfw/src/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../boinc/boinc/samples/gfx_html/browser.cpp \
../boinc/boinc/samples/gfx_html/browserctrl_win.cpp \
../boinc/boinc/samples/gfx_html/browserlog.cpp \
../boinc/boinc/samples/gfx_html/browsermain_win.cpp \
../boinc/boinc/samples/gfx_html/browserwnd_win.cpp \
../boinc/boinc/samples/gfx_html/graphics.cpp \
../boinc/boinc/samples/gfx_html/mongoose.cpp \
../boinc/boinc/samples/gfx_html/vboxwrapper.cpp \
../boinc/boinc/samples/gfx_html/webapi.cpp \
../boinc/boinc/samples/gfx_html/webboincjs.cpp \
../boinc/boinc/samples/gfx_html/webboincpng.cpp \
../boinc/boinc/samples/gfx_html/webindexhtml.cpp \
../boinc/boinc/samples/gfx_html/webserver.cpp 

OBJS += \
./boinc/boinc/samples/gfx_html/browser.o \
./boinc/boinc/samples/gfx_html/browserctrl_win.o \
./boinc/boinc/samples/gfx_html/browserlog.o \
./boinc/boinc/samples/gfx_html/browsermain_win.o \
./boinc/boinc/samples/gfx_html/browserwnd_win.o \
./boinc/boinc/samples/gfx_html/graphics.o \
./boinc/boinc/samples/gfx_html/mongoose.o \
./boinc/boinc/samples/gfx_html/vboxwrapper.o \
./boinc/boinc/samples/gfx_html/webapi.o \
./boinc/boinc/samples/gfx_html/webboincjs.o \
./boinc/boinc/samples/gfx_html/webboincpng.o \
./boinc/boinc/samples/gfx_html/webindexhtml.o \
./boinc/boinc/samples/gfx_html/webserver.o 

CPP_DEPS += \
./boinc/boinc/samples/gfx_html/browser.d \
./boinc/boinc/samples/gfx_html/browserctrl_win.d \
./boinc/boinc/samples/gfx_html/browserlog.d \
./boinc/boinc/samples/gfx_html/browsermain_win.d \
./boinc/boinc/samples/gfx_html/browserwnd_win.d \
./boinc/boinc/samples/gfx_html/graphics.d \
./boinc/boinc/samples/gfx_html/mongoose.d \
./boinc/boinc/samples/gfx_html/vboxwrapper.d \
./boinc/boinc/samples/gfx_html/webapi.d \
./boinc/boinc/samples/gfx_html/webboincjs.d \
./boinc/boinc/samples/gfx_html/webboincpng.d \
./boinc/boinc/samples/gfx_html/webindexhtml.d \
./boinc/boinc/samples/gfx_html/webserver.d 


# Each subdirectory must supply rules for building sources it contributes
boinc/boinc/samples/gfx_html/%.o: ../boinc/boinc/samples/gfx_html/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../lua/src/lapi.c \
../lua/src/lauxlib.c \
../lua/src/lbaselib.c \
../lua/src/lcode.c \
../lua/src/ldblib.c \
../lua/src/ldebug.c \
../lua/src/ldo.c \
../lua/src/ldump.c \
../lua/src/lfunc.c \
../lua/src/lgc.c \
../lua/src/linit.c \
../lua/src/liolib.c \
../lua/src/llex.c \
../lua/src/lmathlib.c \
../lua/src/lmem.c \
../lua/src/loadlib.c \
../lua/src/loadlib_rel.c \
../lua/src/lobject.c \
../lua/src/lopcodes.c \
../lua/src/loslib.c \
../lua/src/lparser.c \
../lua/src/lstate.c \
../lua/src/lstring.c \
../lua/src/lstrlib.c \
../lua/src/ltable.c \
../lua/src/ltablib.c \
../lua/src/ltm.c \
../lua/src/lua.c \
../lua/src/luac.c \
../lua/src/lundump.c \
../lua/src/lvm.c \
../lua/src/lzio.c \
../lua/src/print.c 

OBJS += \
./lua/src/lapi.o \
./lua/src/lauxlib.o \
./lua/src/lbaselib.o \
./lua/src/lcode.o \
./lua/src/ldblib.o \
./lua/src/ldebug.o \
./lua/src/ldo.o \
./lua/src/ldump.o \
./lua/src/lfunc.o \
./lua/src/lgc.o \
./lua/src/linit.o \
./lua/src/liolib.o \
./lua/src/llex.o \
./lua/src/lmathlib.o \
./lua/src/lmem.o \
./lua/src/loadlib.o \
./lua/src/loadlib_rel.o \
./lua/src/lobject.o \
./lua/src/lopcodes.o \
./lua/src/loslib.o \
./lua/src/lparser.o \
./lua/src/lstate.o \
./lua/src/lstring.o \
./lua/src/lstrlib.o \
./lua/src/ltable.o \
./lua/src/ltablib.o \
./lua/src/ltm.o \
./lua/src/lua.o \
./lua/src/luac.o \
./lua/src/lundump.o \
./lua/src/lvm.o \
./lua/src/lzio.o \
./lua/src/print.o 

C_DEPS += \
./lua/src/lapi.d \
./lua/src/lauxlib.d \
./lua/src/lbaselib.d \
./lua/src/lcode.d \
./lua/src/ldblib.d \
./lua/src/ldebug.d \
./lua/src/ldo.d \
./lua/src/ldump.d \
./lua/src/lfunc.d \
./lua/src/lgc.d \
./lua/src/linit.d \
./lua/src/liolib.d \
./lua/src/llex.d \
./lua/src/lmathlib.d \
./lua/src/lmem.d \
./lua/src/loadlib.d \
./lua/src/loadlib_rel.d \
./lua/src/lobject.d \
./lua/src/lopcodes.d \
./lua/src/loslib.d \
./lua/src/lparser.d \
./lua/src/lstate.d \
./lua/src/lstring.d \
./lua/src/lstrlib.d \
./lua/src/ltable.d \
./lua/src/ltablib.d \
./lua/src/ltm.d \
./lua/src/lua.d \
./lua/src/luac.d \
./lua/src/lundump.d \
./lua/src/lvm.d \
./lua/src/lzio.d \
./lua/src/print.d 


# Each subdirectory must supply rules for building sources it contributes
lua/src/%.o: ../lua/src/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



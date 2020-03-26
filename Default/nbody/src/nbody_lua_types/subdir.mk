################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../nbody/src/nbody_lua_types/nbody_lua_body.c \
../nbody/src/nbody_lua_types/nbody_lua_disk.c \
../nbody/src/nbody_lua_types/nbody_lua_dwarf.c \
../nbody/src/nbody_lua_types/nbody_lua_halo.c \
../nbody/src/nbody_lua_types/nbody_lua_histogram_params.c \
../nbody/src/nbody_lua_types/nbody_lua_nbodyctx.c \
../nbody/src/nbody_lua_types/nbody_lua_nbodystate.c \
../nbody/src/nbody_lua_types/nbody_lua_potential.c \
../nbody/src/nbody_lua_types/nbody_lua_spherical.c \
../nbody/src/nbody_lua_types/nbody_lua_type_marshal.c \
../nbody/src/nbody_lua_types/nbody_lua_types.c 

OBJS += \
./nbody/src/nbody_lua_types/nbody_lua_body.o \
./nbody/src/nbody_lua_types/nbody_lua_disk.o \
./nbody/src/nbody_lua_types/nbody_lua_dwarf.o \
./nbody/src/nbody_lua_types/nbody_lua_halo.o \
./nbody/src/nbody_lua_types/nbody_lua_histogram_params.o \
./nbody/src/nbody_lua_types/nbody_lua_nbodyctx.o \
./nbody/src/nbody_lua_types/nbody_lua_nbodystate.o \
./nbody/src/nbody_lua_types/nbody_lua_potential.o \
./nbody/src/nbody_lua_types/nbody_lua_spherical.o \
./nbody/src/nbody_lua_types/nbody_lua_type_marshal.o \
./nbody/src/nbody_lua_types/nbody_lua_types.o 

C_DEPS += \
./nbody/src/nbody_lua_types/nbody_lua_body.d \
./nbody/src/nbody_lua_types/nbody_lua_disk.d \
./nbody/src/nbody_lua_types/nbody_lua_dwarf.d \
./nbody/src/nbody_lua_types/nbody_lua_halo.d \
./nbody/src/nbody_lua_types/nbody_lua_histogram_params.d \
./nbody/src/nbody_lua_types/nbody_lua_nbodyctx.d \
./nbody/src/nbody_lua_types/nbody_lua_nbodystate.d \
./nbody/src/nbody_lua_types/nbody_lua_potential.d \
./nbody/src/nbody_lua_types/nbody_lua_spherical.d \
./nbody/src/nbody_lua_types/nbody_lua_type_marshal.d \
./nbody/src/nbody_lua_types/nbody_lua_types.d 


# Each subdirectory must supply rules for building sources it contributes
nbody/src/nbody_lua_types/%.o: ../nbody/src/nbody_lua_types/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



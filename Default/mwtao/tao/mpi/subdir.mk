################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CXX_SRCS += \
../mwtao/tao/mpi/assign_device.cxx \
../mwtao/tao/mpi/master_worker.cxx \
../mwtao/tao/mpi/mpi_differential_evolution.cxx \
../mwtao/tao/mpi/mpi_genetic_algorithm.cxx \
../mwtao/tao/mpi/mpi_particle_swarm.cxx 

CXX_DEPS += \
./mwtao/tao/mpi/assign_device.d \
./mwtao/tao/mpi/master_worker.d \
./mwtao/tao/mpi/mpi_differential_evolution.d \
./mwtao/tao/mpi/mpi_genetic_algorithm.d \
./mwtao/tao/mpi/mpi_particle_swarm.d 

OBJS += \
./mwtao/tao/mpi/assign_device.o \
./mwtao/tao/mpi/master_worker.o \
./mwtao/tao/mpi/mpi_differential_evolution.o \
./mwtao/tao/mpi/mpi_genetic_algorithm.o \
./mwtao/tao/mpi/mpi_particle_swarm.o 


# Each subdirectory must supply rules for building sources it contributes
mwtao/tao/mpi/%.o: ../mwtao/tao/mpi/%.cxx
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



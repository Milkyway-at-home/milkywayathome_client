################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CXX_SRCS += \
../mwtao/tao/util/hessian.cxx \
../mwtao/tao/util/matrix.cxx \
../mwtao/tao/util/newton_step.cxx \
../mwtao/tao/util/recombination.cxx \
../mwtao/tao/util/regression.cxx \
../mwtao/tao/util/statistics.cxx 

CXX_DEPS += \
./mwtao/tao/util/hessian.d \
./mwtao/tao/util/matrix.d \
./mwtao/tao/util/newton_step.d \
./mwtao/tao/util/recombination.d \
./mwtao/tao/util/regression.d \
./mwtao/tao/util/statistics.d 

OBJS += \
./mwtao/tao/util/hessian.o \
./mwtao/tao/util/matrix.o \
./mwtao/tao/util/newton_step.o \
./mwtao/tao/util/recombination.o \
./mwtao/tao/util/regression.o \
./mwtao/tao/util/statistics.o 


# Each subdirectory must supply rules for building sources it contributes
mwtao/tao/util/%.o: ../mwtao/tao/util/%.cxx
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CXX_SRCS += \
../mwtao/tao/synchronous_algorithms/gradient.cxx \
../mwtao/tao/synchronous_algorithms/line_search.cxx \
../mwtao/tao/synchronous_algorithms/parameter_sweep.cxx \
../mwtao/tao/synchronous_algorithms/synchronous_gradient_descent.cxx \
../mwtao/tao/synchronous_algorithms/synchronous_newton_method.cxx 

CXX_DEPS += \
./mwtao/tao/synchronous_algorithms/gradient.d \
./mwtao/tao/synchronous_algorithms/line_search.d \
./mwtao/tao/synchronous_algorithms/parameter_sweep.d \
./mwtao/tao/synchronous_algorithms/synchronous_gradient_descent.d \
./mwtao/tao/synchronous_algorithms/synchronous_newton_method.d 

OBJS += \
./mwtao/tao/synchronous_algorithms/gradient.o \
./mwtao/tao/synchronous_algorithms/line_search.o \
./mwtao/tao/synchronous_algorithms/parameter_sweep.o \
./mwtao/tao/synchronous_algorithms/synchronous_gradient_descent.o \
./mwtao/tao/synchronous_algorithms/synchronous_newton_method.o 


# Each subdirectory must supply rules for building sources it contributes
mwtao/tao/synchronous_algorithms/%.o: ../mwtao/tao/synchronous_algorithms/%.cxx
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



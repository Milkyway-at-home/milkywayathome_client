################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CXX_SRCS += \
../mwtao/tao/examples/standard_benchmarks.cxx \
../mwtao/tao/examples/standard_benchmarks_db.cxx 

CXX_DEPS += \
./mwtao/tao/examples/standard_benchmarks.d \
./mwtao/tao/examples/standard_benchmarks_db.d 

OBJS += \
./mwtao/tao/examples/standard_benchmarks.o \
./mwtao/tao/examples/standard_benchmarks_db.o 


# Each subdirectory must supply rules for building sources it contributes
mwtao/tao/examples/%.o: ../mwtao/tao/examples/%.cxx
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O2 -g -Wall -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



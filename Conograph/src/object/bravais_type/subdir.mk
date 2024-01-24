
# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
$(ROOT)/bravais_type/BravaisType.cc \
$(ROOT)/bravais_type/BravaisLattice.o \

OBJS += \
object/bravais_type/BravaisType.o \
object/bravais_type/BravaisLattice.o \

DEPS += \
${addprefix object/bravais_type/, \
BravaisType.d \
BravaisLattice.d \
}


# Each subdirectory must supply rules for building sources it contributes
object/bravais_type/%.o: $(ROOT)/bravais_type/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	@echo g++ $(CXXFLAGS) -o$@ $<
	@g++ $(CXXFLAGS) -o$@ $< && \
	echo -n $(@:%.o=%.d) $(dir $@) > $(@:%.o=%.d) && \
	g++ $(PREFLAGS) $(CXXFLAGS) $< >> $(@:%.o=%.d)
	@echo 'Finished building: $<'
	@echo ' '



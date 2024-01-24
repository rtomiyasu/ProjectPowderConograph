
# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
$(ROOT)/LatticeWithSameQ/LatticeMetricTensor.cc \
$(ROOT)/LatticeWithSameQ/LatticeWithSameQ.cc \
$(ROOT)/LatticeWithSameQ/p_out_same_q.cc \

OBJS += \
object/LatticeWithSameQ/LatticeMetricTensor.o \
object/LatticeWithSameQ/LatticeWithSameQ.o \
object/LatticeWithSameQ/p_out_same_q.o \

DEPS += \
${addprefix object/LatticeWithSameQ/, \
LatticeMetricTensor.d \
LatticeWithSameQ.d \
p_out_same_q.d \
}


# Each subdirectory must supply rules for building sources it contributes
object/LatticeWithSameQ/%.o: $(ROOT)/LatticeWithSameQ/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	@echo g++ $(CXXFLAGS) -o$@ $<
	@g++ $(CXXFLAGS) -o$@ $< && \
	echo -n $(@:%.o=%.d) $(dir $@) > $(@:%.o=%.d) && \
	g++ $(PREFLAGS) $(CXXFLAGS) $< >> $(@:%.o=%.d)
	@echo 'Finished building: $<'
	@echo ' '



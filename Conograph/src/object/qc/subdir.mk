
# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
$(ROOT)/qc/gather_qcal2.cc \
$(ROOT)/qc/reflection_conditions.cc \
$(ROOT)/qc/p_out_space_group_dtm.cc \

OBJS += \
object/qc/gather_qcal2.o \
object/qc/reflection_conditions.o \
object/qc/p_out_space_group_dtm.o \

DEPS += \
${addprefix object/qc/, \
gather_qcal2.d \
reflection_conditions.d \
p_out_space_group_dtm.d \
}


# Each subdirectory must supply rules for building sources it contributes
object/qc/%.o: $(ROOT)/qc/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	@echo g++ $(CXXFLAGS) -o$@ $<
	@g++ $(CXXFLAGS) -o$@ $< && \
	echo -n $(@:%.o=%.d) $(dir $@) > $(@:%.o=%.d) && \
	g++ $(PREFLAGS) $(CXXFLAGS) $< >> $(@:%.o=%.d)
	@echo 'Finished building: $<'
	@echo ' '



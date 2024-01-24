
# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
$(ROOT)/utility_rw_param/I_ReadData.cc \
$(ROOT)/utility_rw_param/RWParam_void.cc \

OBJS += \
object/utility_rw_param/I_ReadData.o \
object/utility_rw_param/RWParam_void.o \

DEPS += \
${addprefix object/utility_rw_param/, \
I_ReadData.d \
RWParam_void.d \
}


# Each subdirectory must supply rules for building sources it contributes
object/utility_rw_param/%.o: $(ROOT)/utility_rw_param/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	@echo g++ $(CXXFLAGS) -o$@ $<
	@g++ $(CXXFLAGS) -o$@ $< && \
	echo -n $(@:%.o=%.d) $(dir $@) > $(@:%.o=%.d) && \
	g++ $(PREFLAGS) $(CXXFLAGS) $< >> $(@:%.o=%.d)
	@echo 'Finished building: $<'
	@echo ' '



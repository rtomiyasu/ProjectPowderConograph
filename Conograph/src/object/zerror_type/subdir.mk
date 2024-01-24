# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
$(ROOT)/zerror_type/error_mes.cc \

OBJS += \
object/zerror_type/error_mes.o \

DEPS += \
${addprefix object/zerror_type/, \
error_mes.d \
}


# Each subdirectory must supply rules for building sources it contributes
object/zerror_type/%.o: $(ROOT)/zerror_type/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	@echo g++ $(CXXFLAGS) -c -fmessage-length=0 -o$@ $<
	@g++ $(CXXFLAGS) -c -fmessage-length=0 -o$@ $< && \
	echo -n $(@:%.o=%.d) $(dir $@) > $(@:%.o=%.d) && \
	g++ -MM -MG -P -w $(CXXFLAGS) -c -fmessage-length=0  $< >> $(@:%.o=%.d)
	@echo 'Finished building: $<'
	@echo ' '



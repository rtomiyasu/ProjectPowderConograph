# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
$(ROOT)/zparam/etype_ID.cc \

OBJS += \
object/zparam/etype_ID.o \

DEPS += \
${addprefix object/zparam/, \
etype_ID.d \
}


# Each subdirectory must supply rules for building sources it contributes
object/zparam/%.o: $(ROOT)/zparam/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	@echo g++ $(CXXFLAGS) -c -fmessage-length=0 -o$@ $<
	@g++ $(CXXFLAGS) -c -fmessage-length=0 -o$@ $< && \
	echo -n $(@:%.o=%.d) $(dir $@) > $(@:%.o=%.d) && \
	g++ -MM -MG -P -w $(CXXFLAGS) -c -fmessage-length=0  $< >> $(@:%.o=%.d)
	@echo 'Finished building: $<'
	@echo ' '



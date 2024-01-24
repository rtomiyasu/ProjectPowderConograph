
# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
$(ROOT)/centring_type/CentringType.cc \
$(ROOT)/centring_type/bravais_lat.cc \

OBJS += \
object/centring_type/CentringType.o \
object/centring_type/bravais_lat.o \

DEPS += \
${addprefix object/centring_type/, \
CentringType.d \
bravais_lat.d \
}


# Each subdirectory must supply rules for building sources it contributes
object/centring_type/%.o: $(ROOT)/centring_type/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	@echo g++ $(CXXFLAGS) -o$@ $<
	@g++ $(CXXFLAGS) -o$@ $< && \
	echo -n $(@:%.o=%.d) $(dir $@) > $(@:%.o=%.d) && \
	g++ $(PREFLAGS) $(CXXFLAGS) $< >> $(@:%.o=%.d)
	@echo 'Finished building: $<'
	@echo ' '



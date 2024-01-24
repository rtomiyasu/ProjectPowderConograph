
# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
$(ROOT)/laue_group/LaueGroup.cc \
$(ROOT)/laue_group/laue_gp.cc \

OBJS += \
object/laue_group/LaueGroup.o \
object/laue_group/laue_gp.o \

DEPS += \
${addprefix object/laue_group/, \
LaueGroup.d \
laue_gp.d \
}


# Each subdirectory must supply rules for building sources it contributes
object/laue_group/%.o: $(ROOT)/laue_group/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	@echo g++ $(CXXFLAGS) -o$@ $<
	@g++ $(CXXFLAGS) -o$@ $< && \
	echo -n $(@:%.o=%.d) $(dir $@) > $(@:%.o=%.d) && \
	g++ $(PREFLAGS) $(CXXFLAGS) $< >> $(@:%.o=%.d)
	@echo 'Finished building: $<'
	@echo ' '



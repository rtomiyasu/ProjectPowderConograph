# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
$(ROOT)/levenberg_marquardt/LemarqMethod.cc \
$(ROOT)/levenberg_marquardt/LemarqMethod_Marquardt_Original.cc \
$(ROOT)/levenberg_marquardt/SVdcmp.cc \

OBJS += \
object/levenberg_marquardt/LemarqMethod.o \
object/levenberg_marquardt/LemarqMethod_Marquardt_Original.o \
object/levenberg_marquardt/SVdcmp.o \

DEPS += \
${addprefix object/levenberg_marquardt/, \
LemarqMethod.d \
LemarqMethod_Marquardt_Original.d \
SVdcmp.d \
}


# Each subdirectory must supply rules for building sources it contributes
object/levenberg_marquardt/%.o: $(ROOT)/levenberg_marquardt/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	@echo g++ $(CXXFLAGS) -o$@ $<
	@g++ $(CXXFLAGS) -o$@ $< && \
	echo -n $(@:%.o=%.d) $(dir $@) > $(@:%.o=%.d) && \
	g++ $(PREFLAGS) $(CXXFLAGS) $< >> $(@:%.o=%.d)
	@echo 'Finished building: $<'
	@echo ' '



# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
$(ROOT)/zlog/rlog.cc \
$(ROOT)/zlog/zlog.cc \

OBJS += \
object/zlog/rlog.o \
object/zlog/zlog.o \

DEPS += \
${addprefix object/zlog/, \
rlog.d \
zlog.d \
}


# Each subdirectory must supply rules for building sources it contributes
object/zlog/%.o: $(ROOT)/zlog/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	@echo g++ $(CXXFLAGS) -c -fmessage-length=0 -o$@ $<
	@g++ $(CXXFLAGS) -c -fmessage-length=0 -o$@ $< && \
	echo -n $(@:%.o=%.d) $(dir $@) > $(@:%.o=%.d) && \
	g++ -MM -MG -P -w $(CXXFLAGS) -c -fmessage-length=0  $< >> $(@:%.o=%.d)
	@echo 'Finished building: $<'
	@echo ' '

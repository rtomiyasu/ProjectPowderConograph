
# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
$(ROOT)/model_function/profile_function/global_function/GlbBraggDiffract.cc \
$(ROOT)/model_function/profile_function/global_function/GlbPolynomialConv.cc \
$(ROOT)/model_function/profile_function/global_function/IGlobalFunc.cc \
$(ROOT)/model_function/profile_function/global_function/PeakShiftFunc.cc \

OBJS += \
object/model_function/profile_function/global_function/GlbBraggDiffract.o \
object/model_function/profile_function/global_function/GlbPolynomialConv.o \
object/model_function/profile_function/global_function/IGlobalFunc.o \
object/model_function/profile_function/global_function/PeakShiftFunc.o \

DEPS += \
${addprefix object/model_function/profile_function/global_function/, \
GlbBraggDiffract.d \
GlbPolynomialConv.d \
IGlobalFunc.d \
PeakShiftFunc.d \
}


# Each subdirectory must supply rules for building sources it contributes
object/model_function/profile_function/global_function/%.o: $(ROOT)/model_function/profile_function/global_function/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	@echo g++ $(CXXFLAGS) -c -fmessage-length=0 -o$@ $<
	@g++ $(CXXFLAGS) -c -fmessage-length=0 -o$@ $< && \
	echo -n $(@:%.o=%.d) $(dir $@) > $(@:%.o=%.d) && \
	g++ -MM -MG -P -w $(CXXFLAGS) -c -fmessage-length=0  $< >> $(@:%.o=%.d)
	@echo 'Finished building: $<'
	@echo ' '



# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
$(ROOT)/model_function/LatticeDistanceModel.cc \
$(ROOT)/model_function/MarquardtFmodelBase.cc \
$(ROOT)/model_function/PeakPosModel.cc \

OBJS += \
object/model_function/LatticeDistanceModel.o \
object/model_function/MarquardtFmodelBase.o \
object/model_function/PeakPosModel.o \

DEPS += \
${addprefix object/model_function/, \
LatticeDistanceModel.d \
MarquardtFmodelBase.d \
PeakPosModel.d \
}


# Each subdirectory must supply rules for building sources it contributes
object/model_function/%.o: $(ROOT)/model_function/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	@echo g++ $(CXXFLAGS) -o$@ $<
	@g++ $(CXXFLAGS) -o$@ $< && \
	echo -n $(@:%.o=%.d) $(dir $@) > $(@:%.o=%.d) && \
	g++ $(PREFLAGS) $(CXXFLAGS) $< >> $(@:%.o=%.d)
	@echo 'Finished building: $<'
	@echo ' '



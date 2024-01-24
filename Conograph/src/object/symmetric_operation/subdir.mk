
# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
$(ROOT)/symmetric_operation/MillerIndex.cc \
$(ROOT)/symmetric_operation/S1.cc \
$(ROOT)/symmetric_operation/StringS1.cc \
$(ROOT)/symmetric_operation/SymmetricOperation.cc \
$(ROOT)/symmetric_operation/translation_vector.cc \

OBJS += \
object/symmetric_operation/MillerIndex.o \
object/symmetric_operation/S1.o \
object/symmetric_operation/StringS1.o \
object/symmetric_operation/SymmetricOperation.o \
object/symmetric_operation/translation_vector.o \

DEPS += \
${addprefix object/symmetric_operation/, \
MillerIndex.d \
S1.d \
StringS1.d \
SymmetricOperation.d \
translation_vector.d \
}


# Each subdirectory must supply rules for building sources it contributes
object/symmetric_operation/%.o: $(ROOT)/symmetric_operation/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	@echo g++ $(CXXFLAGS) -o$@ $<
	@g++ $(CXXFLAGS) -o$@ $< && \
	echo -n $(@:%.o=%.d) $(dir $@) > $(@:%.o=%.d) && \
	g++ $(PREFLAGS) $(CXXFLAGS) $< >> $(@:%.o=%.d)
	@echo 'Finished building: $<'
	@echo ' '



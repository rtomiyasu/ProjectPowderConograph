
# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
$(ROOT)/utility_func/covar_matrix.cc \
$(ROOT)/utility_func/gcd.cc \
$(ROOT)/utility_func/lattice_constant.cc \
$(ROOT)/utility_func/zmath.cc \
$(ROOT)/utility_func/stopx.cc \
$(ROOT)/utility_func/zstring.cc \

OBJS += \
object/utility_func/covar_matrix.o \
object/utility_func/gcd.o \
object/utility_func/lattice_constant.o \
object/utility_func/zmath.o \
object/utility_func/stopx.o \
object/utility_func/zstring.o \

DEPS += \
${addprefix object/utility_func/, \
covar_matrix.d \
gcd.d \
lattice_constant.d \
zmath.d \
stopx.d \
zstring.d \
}


# Each subdirectory must supply rules for building sources it contributes
object/utility_func/%.o: $(ROOT)/utility_func/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	@echo g++ $(CXXFLAGS) -o$@ $<
	@g++ $(CXXFLAGS) -o$@ $< && \
	echo -n $(@:%.o=%.d) $(dir $@) > $(@:%.o=%.d) && \
	g++ $(PREFLAGS) $(CXXFLAGS) $< >> $(@:%.o=%.d)
	@echo 'Finished building: $<'
	@echo ' '



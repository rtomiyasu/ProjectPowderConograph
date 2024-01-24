
# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
$(ROOT)/point_group/coset_representative_data.cc \
$(ROOT)/point_group/PGNormalSeriesTray.cc \
$(ROOT)/point_group/point_gp_data.cc \

OBJS += \
object/point_group/coset_representative_data.o \
object/point_group/PGNormalSeriesTray.o \
object/point_group/point_gp_data.o \

DEPS += \
${addprefix object/point_group/, \
coset_representative_data.d \
PGNormalSeriesTray.d \
point_gp_data.d \
}


# Each subdirectory must supply rules for building sources it contributes
object/point_group/%.o: $(ROOT)/point_group/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	@echo g++ $(CXXFLAGS) -o$@ $<
	@g++ $(CXXFLAGS) -o$@ $< && \
	echo -n $(@:%.o=%.d) $(dir $@) > $(@:%.o=%.d) && \
	g++ $(PREFLAGS) $(CXXFLAGS) $< >> $(@:%.o=%.d)
	@echo 'Finished building: $<'
	@echo ' '



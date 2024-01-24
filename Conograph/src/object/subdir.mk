
# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
$(ROOT)/chToqValue.cc \
$(ROOT)/ControlFile.cc \
$(ROOT)/ControlParam.cc \
$(ROOT)/p_out_indexing.cc \
$(ROOT)/GenerateRelation.cc \
$(ROOT)/IndexingLattice.cc \
$(ROOT)/PeakPosData.cc \
$(ROOT)/indexing_func_dim2.cc \
$(ROOT)/indexing_func_dim3.cc \
$(ROOT)/main_indexing.cc \
$(ROOT)/SortingLattice.cc \

OBJS += \
object/chToqValue.o \
object/ControlFile.o \
object/ControlParam.o \
object/p_out_indexing.o \
object/GenerateRelation.o \
object/IndexingLattice.o \
object/PeakPosData.o \
object/indexing_func_dim2.o \
object/indexing_func_dim3.o \
object/main_indexing.o \
object/SortingLattice.o \

DEPS += \
${addprefix object/, \
chToqValue.d \
ControlFile.d \
ControlParam.d \
p_out_indexing.d \
GenerateRelation.d \
IndexingLattice.d \
PeakPosData.d \
indexing_func_dim2.d \
indexing_func_dim3.d \
main_indexing.d \
SortingLattice.d \
}


# Each subdirectory must supply rules for building sources it contributes
object/%.o: $(ROOT)/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	@echo g++ $(CXXFLAGS) -o$@ $<
	@g++ $(CXXFLAGS) -o$@ $< && \
	echo -n $(@:%.o=%.d) $(dir $@) > $(@:%.o=%.d) && \
	g++ $(PREFLAGS) $(CXXFLAGS) $< >> $(@:%.o=%.d)
	@echo 'Finished building: $<'
	@echo ' '



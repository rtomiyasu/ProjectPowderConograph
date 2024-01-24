
# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
$(ROOT)/lattice_symmetry/gather_q_of_3D_lattice.cc \
$(ROOT)/lattice_symmetry/gather_q_of_Ndim_lattice.cc \
$(ROOT)/lattice_symmetry/LatticeFigureOfMerit.cc \
$(ROOT)/lattice_symmetry/LatticeFigureOfMeritZeroShift.cc \
$(ROOT)/lattice_symmetry/LatticeFigureOfMeritToCheckSymmetry.cc \
$(ROOT)/lattice_symmetry/LatticeFigureOfMeritToDisplay.cc \
$(ROOT)/lattice_symmetry/OutputInfo.cc \
$(ROOT)/lattice_symmetry/ReducedLatticeToCheckBravais.cc \
$(ROOT)/lattice_symmetry/ReducedVCLatticeToCheckBravais.cc \
$(ROOT)/lattice_symmetry/ReducedLatticeToCheckEquiv.cc \
$(ROOT)/lattice_symmetry/VCLatticeFigureOfMeritToCheckSymmetry.cc \

OBJS += \
object/lattice_symmetry/gather_q_of_3D_lattice.o \
object/lattice_symmetry/gather_q_of_Ndim_lattice.o \
object/lattice_symmetry/LatticeFigureOfMerit.o \
object/lattice_symmetry/LatticeFigureOfMeritZeroShift.o \
object/lattice_symmetry/LatticeFigureOfMeritToCheckSymmetry.o \
object/lattice_symmetry/LatticeFigureOfMeritToDisplay.o \
object/lattice_symmetry/OutputInfo.o \
object/lattice_symmetry/ReducedLatticeToCheckBravais.o \
object/lattice_symmetry/ReducedVCLatticeToCheckBravais.o \
object/lattice_symmetry/ReducedLatticeToCheckEquiv.o \
object/lattice_symmetry/VCLatticeFigureOfMeritToCheckSymmetry.o \

DEPS += \
${addprefix object/lattice_symmetry/, \
gather_q_of_3D_lattice.d \
gather_q_of_Ndim_lattice.d \
LatticeFigureOfMerit.d \
LatticeFigureOfMeritZeroShift.d \
LatticeFigureOfMeritToCheckSymmetry.d \
LatticeFigureOfMeritToDisplay.d \
lattice_symmetry/OutputInfo.d \
ReducedLatticeToCheckBravais.d \
ReducedVCLatticeToCheckBravais.d \
ReducedLatticeToCheckEquiv.d \
VCLatticeFigureOfMeritToCheckSymmetry.d \
}


# Each subdirectory must supply rules for building sources it contributes
object/lattice_symmetry/%.o: $(ROOT)/lattice_symmetry/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	@echo g++ $(CXXFLAGS) -o$@ $<
	@g++ $(CXXFLAGS) -o$@ $< && \
	echo -n $(@:%.o=%.d) $(dir $@) > $(@:%.o=%.d) && \
	g++ $(PREFLAGS) $(CXXFLAGS) $< >> $(@:%.o=%.d)
	@echo 'Finished building: $<'
	@echo ' '



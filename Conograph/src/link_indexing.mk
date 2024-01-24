# object files linked for building exe files

# ROOT
OBJS_INDEXING += \
object/ControlFile.o \
object/ControlParam.o \
object/chToqValue.o \
object/p_out_indexing.o \
object/GenerateRelation.o \
object/IndexingLattice.o \
object/PeakPosData.o \
object/indexing_func_dim2.o \
object/indexing_func_dim3.o \
object/main_indexing.o \
object/SortingLattice.o \

OBJS_INDEXING += \
object/zerror_type/error_mes.o \

OBJS_INDEXING += \
object/utility_func/covar_matrix.o \
object/utility_func/gcd.o \
object/utility_func/lattice_constant.o \
object/utility_func/zmath.o \
object/utility_func/zstring.o \
object/utility_func/stopx.o \


OBJS_INDEXING += \
object/utility_rw_param/I_ReadData.o \
object/utility_rw_param/RWParam_void.o \


OBJS_INDEXING += \
object/utility_data_structure/Bud.o \
object/utility_data_structure/Bud2.o \
object/utility_data_structure/Node3.o \
object/utility_data_structure/NodeB.o \
object/utility_data_structure/SymMatWCovar.o \
object/utility_data_structure/TreeLattice.o \
object/utility_data_structure/VCData.o \


# bravais_lattice
OBJS_INDEXING += \
object/centring_type/CentringType.o \
object/centring_type/bravais_lat.o \


# bravais_type
OBJS_INDEXING += \
object/bravais_type/BravaisType.o \
object/bravais_type/BravaisLattice.o \


#laue_group
OBJS_INDEXING += \
object/laue_group/LaueGroup.o \
object/laue_group/laue_gp.o \


# point_group
OBJS_INDEXING += \
object/point_group/coset_representative_data.o \
object/point_group/PGNormalSeriesTray.o \
object/point_group/point_gp_data.o \


# symmetric_operation
OBJS_INDEXING += \
object/symmetric_operation/MillerIndex.o \
object/symmetric_operation/S1.o \
object/symmetric_operation/SymmetricOperation.o \
object/symmetric_operation/translation_vector.o \
object/symmetric_operation/StringS1.o \


# lattice_symmetry
OBJS_INDEXING += \
object/lattice_symmetry/gather_q_of_3D_lattice.o \
object/lattice_symmetry/gather_q_of_Ndim_lattice.o \
object/lattice_symmetry/LatticeFigureOfMerit.o \
object/lattice_symmetry/LatticeFigureOfMeritToCheckSymmetry.o \
object/lattice_symmetry/LatticeFigureOfMeritToDisplay.o \
object/lattice_symmetry/LatticeFigureOfMeritZeroShift.o \
object/lattice_symmetry/OutputInfo.o \
object/lattice_symmetry/ReducedLatticeToCheckBravais.o \
object/lattice_symmetry/ReducedVCLatticeToCheckBravais.o \
object/lattice_symmetry/ReducedLatticeToCheckEquiv.o \
object/lattice_symmetry/VCLatticeFigureOfMeritToCheckSymmetry.o \

OBJS_INDEXING += \
object/LatticeWithSameQ/LatticeMetricTensor.o \
object/LatticeWithSameQ/LatticeWithSameQ.o \
object/LatticeWithSameQ/p_out_same_q.o \


OBJS_INDEXING += \
object/levenberg_marquardt/LemarqMethod.o \
object/levenberg_marquardt/LemarqMethod_Marquardt_Original.o \
object/levenberg_marquardt/SVdcmp.o \

OBJS_INDEXING += \
object/model_function/LatticeDistanceModel.o \
object/model_function/MarquardtFmodelBase.o \
object/model_function/PeakPosModel.o \

# model_function/profile_function/global_function
OBJS_INDEXING += \
object/model_function/profile_function/global_function/GlbBraggDiffract.o \
object/model_function/profile_function/global_function/GlbPolynomialConv.o \
object/model_function/profile_function/global_function/IGlobalFunc.o \
object/model_function/profile_function/global_function/PeakShiftFunc.o \

# qc
OBJS_INDEXING += \
object/qc/reflection_conditions.o \
object/qc/gather_qcal2.o \
object/qc/p_out_space_group_dtm.o \


# zlog
OBJS_INDEXING += \
object/zlog/rlog.o \
object/zlog/zlog.o \

# zparam
OBJS_INDEXING += \
object/zparam/etype_ID.o \

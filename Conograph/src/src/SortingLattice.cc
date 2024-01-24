/*
 * The MIT License

   Conograph (powder auto-indexing program)

Copyright (c) <2012> <Ryoko Oishi-Tomiyasu, KEK>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
 *
 */
#ifdef _OPENMP
 # include <omp.h>
#endif
#include "utility_func/chToDouble.hh"
#include "utility_func/transform_sym_matrix.hh"
#include "utility_lattice_reduction/put_Selling_reduced_lattice.hh"
#include "lattice_symmetry/VCLatticeFigureOfMeritToCheckSymmetry.hh"
#include "lattice_symmetry/ReducedLatticeToCheckEquiv.hh"
#include "zerror_type/error_out.hh"
#include "zlog/zlog.hh"
#include "ControlParam.hh"
#include "utility_func/stopx.hh"
#include "SortingLattice.hh"

const bool SortingLattice::m_DoesPrudentSymSearch = false;
const Double SortingLattice::m_cv2 = 0.5;

SortingLattice::SortingLattice()
{
	for(Int4 i=0; i<NUM_LS; i++)
	{
		OutputSymmetry[i] = false;
		JudgeSymmetry[i] = false;
	}
	
   	m_resol = 0.0;
	m_num_ref_figure_of_merit = 20;
	m_etype_peak_shift = kPeakShiftFunction_Type0;
	m_WlengthX = 1.54056;
}


SortingLattice::~SortingLattice()
{
}


// Set the member variables.
void SortingLattice::setParam(const ControlParam& cont)
{
	OutputSymmetry[(size_t)Triclinic] = cont.putOutputSymmetry(Triclinic);
	JudgeSymmetry[(size_t)Triclinic] = false;
	for(Int4 i=1; i<NUM_LS; i++)
	{
		OutputSymmetry[i] = cont.putOutputSymmetry(eBravaisType(i));
		JudgeSymmetry[i] = cont.putOutputSymmetry(eBravaisType(i));
	}

	if( JudgeSymmetry[(size_t)Cubic_P] )
	{
		JudgeSymmetry[(size_t)Tetragonal_P] = true;
	}
	if( JudgeSymmetry[(size_t)Hexagonal] )
	{
		JudgeSymmetry[(size_t)Monoclinic_P] = true;
	}
	if( JudgeSymmetry[(size_t)Tetragonal_P] )
	{
		JudgeSymmetry[(size_t)Orthorhombic_P] = true;
	}
	if( JudgeSymmetry[(size_t)Orthorhombic_P] )
	{
		JudgeSymmetry[(size_t)Monoclinic_P] = true;
	}
	
	if( JudgeSymmetry[(size_t)Orthorhombic_C] )
	{
		JudgeSymmetry[(size_t)Monoclinic_B] = true;
	}

	if( JudgeSymmetry[(size_t)Cubic_I] )
	{
		JudgeSymmetry[(size_t)Tetragonal_I] = true;
	}
	if( JudgeSymmetry[(size_t)Tetragonal_I] )
	{
		JudgeSymmetry[(size_t)Orthorhombic_I] = true;
	}

	if( JudgeSymmetry[(size_t)Cubic_F] )
	{
		JudgeSymmetry[(size_t)Orthorhombic_F] = true;
	}

	m_resol = cont.putResolution();
	m_num_ref_figure_of_merit = cont.putNumberOfReflectionsForFigureOfMerit();
	m_etype_peak_shift = cont.putPeakShiftFunctionType();
	m_WlengthX = cont.putWaveLength();

	const vector<Double>& peak_shift_param_rad = cont.putPeakShiftParamRadian();
	const Int4 param_num = peak_shift_param_rad.size();
	assert( m_etype_peak_shift != kPeakShiftFunction_Type0 || param_num == 0 );
	assert( m_etype_peak_shift != kPeakShiftFunction_Type1 || param_num == 1 );

	m_peak_shift_param_rad.resize(param_num);
	for(Int4 i=0; i<param_num; i++)
	{
		m_peak_shift_param_rad[i] = peak_shift_param_rad[i];
	}
}




void SortingLattice::putCentringTypes(const ReducedVCLatticeToCheckBravais& RLCB,
		const VCLatticeFigureOfMeritToCheckSymmetry& lattice_original,
		const BravaisType& brat,
		vector<VCLatticeFigureOfMeritToCheckSymmetry>& lattice_result) const
{
	lattice_result.clear();

	const map< SymMat<VCData>, NRMat<Int4> >& S_red_tray = RLCB.checkCentringType(brat);
	if( S_red_tray.empty() ) return;

	// The lattice of RLCB has at least the symmetry given by eblat.
	SymMat<VCData> S_super(4);
	NRMat<Int4> trans_mat(4,3);

	for(map< SymMat<VCData>, NRMat<Int4> >::const_iterator it=S_red_tray.begin(); it!=S_red_tray.end(); it++)
	{
		S_super = transform_sym_matrix(it->second, it->first);
		trans_mat = identity_matrix<Int4>(4);

		// S_super = trans_mat * it->second * it->first * Transpose(trans_mat * it->second) is Delone reduced.
		if( !put_Selling_reduced_dim_less_than_4(S_super, trans_mat) )
		{
			assert( false );
		}
		moveSmallerDiagonalLeftUpper(S_super, trans_mat);

		lattice_result.push_back( VCLatticeFigureOfMeritToCheckSymmetry( brat, SymMat43_VCData(it->first, mprod(trans_mat, it->second) ),
									lattice_original.putLatticeFigureOfMerit().putPeakShiftFunctionType(),
									lattice_original.putLatticeFigureOfMerit().putWaveLength(),
									lattice_original.putLatticeFigureOfMerit().putPeakShiftParamRadian() ) );
	}
}






void SortingLattice::putLatticeCandidatesForTriclinic(const vector<SymMat43_VCData>& S_super,
		const Double& MIN_NormM,
		const Double& MIN_RevM,
		vector<VCLatticeFigureOfMeritToCheckSymmetry>& lattice_result_tri) const
{
	const Int4 num_topo = S_super.size();
	lattice_result_tri.clear();

/*   2011.10.19 VIC Tamura */
Int4 LOOP_COUNTER = 0;
	
#ifdef _OPENMP
	#pragma omp parallel
#endif
	{
		vector< VecDat3<Int4> > closest_hkl_tray;
		vector<bool> is_cal_Q_observed_tray;
		vector<VCLatticeFigureOfMeritToCheckSymmetry> latFOM_tray;

#ifdef _OPENMP
		#pragma omp for
#endif
		for(Int4 n=0; n<num_topo; n++)
	   	{
/*   2011.10.19 VIC Tamura */
SET_PROGRESS(100*(LOOP_COUNTER++)/num_topo, 65, 1); // critical, but works
if(IS_CANSELED()) continue;

			VCLatticeFigureOfMeritToCheckSymmetry latFOM(BravaisType( pair<eCentringType, ePointGroup>(Prim, Ci) ), S_super[n],
											m_etype_peak_shift, m_WlengthX, m_peak_shift_param_rad);

			latFOM.setFigureOfMerit(m_num_ref_figure_of_merit,
										VCData::putPeakQData(),
										closest_hkl_tray, is_cal_Q_observed_tray);

			LatticeFigureOfMeritZeroShift latFOM2 = latFOM.putLatticeFigureOfMerit();
			pair<bool, ZErrorMessage> ans = latFOM2.fitLatticeParameterLinear(VCData::putPeakQData(),
												closest_hkl_tray, is_cal_Q_observed_tray, false);

			if( ans.first )
			{
				assert( latFOM.putLatticeFigureOfMerit().putFiguresOfMerit().putNumberOfReflectionsForFigureOfMerit() > 0 );
				latFOM2.setFigureOfMerit(latFOM.putLatticeFigureOfMerit().putFiguresOfMerit().putNumberOfReflectionsForFigureOfMerit(),
											VCData::putPeakQData(),
											closest_hkl_tray, is_cal_Q_observed_tray);
				if( LatticeFigureOfMerit::cmpFOMWu( latFOM2, latFOM.putLatticeFigureOfMerit() ) )
				{
					latFOM.setLatticeFigureOfMerit(latFOM2);
				}
			}
			const LatticeFigureOfMerit::SetOfFigureOfMerit& setFOM = latFOM.putLatticeFigureOfMerit().putFiguresOfMerit();
			if( setFOM.putFigureOfMeritWu() < MIN_NormM ) continue;
			if( setFOM.putReversedFigureOfMerit() < MIN_RevM ) continue;

#ifdef _OPENMP
		#pragma omp critical
#endif
			{
				lattice_result_tri.push_back( latFOM );
			}
	   	}
	}

/*   2011.10.19 VIC Tamura */
CHECK_INTERRUPTION();

//	sort( lattice_result_tri.begin(), lattice_result_tri.end() );
}


void SortingLattice::putLatticeCandidatesForEachBravaisTypes(
					const Double& MIN_NormM,
					const Double& MIN_RevM,
					const Int4& MAX_SIZE,
					const eABCaxis& abc_axis,
					const eRHaxis& rh_axis,
					vector<VCLatticeFigureOfMeritToCheckSymmetry> lattice_result[NUM_LS]) const
{
	try{

	for(Int4 i=1; i<NUM_LS; i++)
	{
		lattice_result[i].clear();
	}
	vector<VCLatticeFigureOfMeritToCheckSymmetry>& lattice_result_tri = lattice_result[(size_t)Triclinic];
	vector<VCLatticeFigureOfMeritToCheckSymmetry>& lattice_result_mono_P = lattice_result[(size_t)Monoclinic_P];
	vector<VCLatticeFigureOfMeritToCheckSymmetry>& lattice_result_mono_B = lattice_result[(size_t)Monoclinic_B];
	vector<VCLatticeFigureOfMeritToCheckSymmetry>& lattice_result_ortho_P = lattice_result[(size_t)Orthorhombic_P];
	vector<VCLatticeFigureOfMeritToCheckSymmetry>& lattice_result_ortho_B = lattice_result[(size_t)Orthorhombic_C];
	vector<VCLatticeFigureOfMeritToCheckSymmetry>& lattice_result_ortho_I = lattice_result[(size_t)Orthorhombic_I];
	vector<VCLatticeFigureOfMeritToCheckSymmetry>& lattice_result_ortho_F = lattice_result[(size_t)Orthorhombic_F];
	vector<VCLatticeFigureOfMeritToCheckSymmetry>& lattice_result_tetra_P = lattice_result[(size_t)Tetragonal_P];
	vector<VCLatticeFigureOfMeritToCheckSymmetry>& lattice_result_tetra_I = lattice_result[(size_t)Tetragonal_I];
	vector<VCLatticeFigureOfMeritToCheckSymmetry>& lattice_result_rhom = lattice_result[(size_t)Rhombohedral];
	vector<VCLatticeFigureOfMeritToCheckSymmetry>& lattice_result_hex = lattice_result[(size_t)Hexagonal];
	vector<VCLatticeFigureOfMeritToCheckSymmetry>& lattice_result_cubic_P = lattice_result[(size_t)Cubic_P];
	vector<VCLatticeFigureOfMeritToCheckSymmetry>& lattice_result_cubic_I = lattice_result[(size_t)Cubic_I];
	vector<VCLatticeFigureOfMeritToCheckSymmetry>& lattice_result_cubic_F = lattice_result[(size_t)Cubic_F];

	const Int4 num_tri = lattice_result_tri.size();

/*   2011.10.19 VIC Tamura */
Int4 LOOP_COUNTER = 0;

#ifdef _OPENMP
	#pragma omp parallel
#endif
	{
		vector<VCLatticeFigureOfMeritToCheckSymmetry> latFOM_tray;

#ifdef _OPENMP
		#pragma omp for
#endif
		for(Int4 n=0; n<num_tri; n++)
	   	{
/*   2011.10.19 VIC Tamura */
SET_PROGRESS(99*(LOOP_COUNTER++)/num_tri, 66, 30); // critical, but works
if(IS_CANSELED()) continue;

			VCLatticeFigureOfMeritToCheckSymmetry& latFOM = lattice_result_tri[n];
			latFOM.setLabel(n+1);

			const ReducedVCLatticeToCheckBravais RLCB(abc_axis, rh_axis, m_DoesPrudentSymSearch, m_cv2, latFOM.putInitialSellingReducedForm() );

			if( JudgeSymmetry[Monoclinic_B] )
			{
				putCentringTypes(RLCB, latFOM, BravaisType( put_monoclinic_b_type(abc_axis) ), latFOM_tray);
#ifdef _OPENMP
				#pragma omp critical(monoB)
#endif
				{
					lattice_result_mono_B.insert(lattice_result_mono_B.end(), latFOM_tray.begin(), latFOM_tray.end());
				}
			}
			if( JudgeSymmetry[Orthorhombic_I] )
			{
				putCentringTypes(RLCB, latFOM, BravaisType( pair<eCentringType, ePointGroup>(Inner, D2h) ), latFOM_tray);
#ifdef _OPENMP
				#pragma omp critical(orthoI)
#endif
				{
					lattice_result_ortho_I.insert(lattice_result_ortho_I.end(), latFOM_tray.begin(), latFOM_tray.end());
				}
			}
			if( JudgeSymmetry[Orthorhombic_F] )
			{
				putCentringTypes(RLCB, latFOM, BravaisType( pair<eCentringType, ePointGroup>(Face, D2h) ), latFOM_tray);
#ifdef _OPENMP
				#pragma omp critical(orthoF)
#endif
				{
					lattice_result_ortho_F.insert(lattice_result_ortho_F.end(), latFOM_tray.begin(), latFOM_tray.end());
				}
			}
			if( JudgeSymmetry[Rhombohedral] )
			{
				putCentringTypes(RLCB, latFOM, BravaisType( put_rhombohedral_type(rh_axis) ), latFOM_tray);
#ifdef _OPENMP
				#pragma omp critical(rhom)
#endif
				{
					lattice_result_rhom.insert(lattice_result_rhom.end(), latFOM_tray.begin(), latFOM_tray.end());
				}
			}

			if( JudgeSymmetry[Monoclinic_P] )
			{
				latFOM.putLatticesOfHigherSymmetry(put_monoclinic_p_type(abc_axis), m_cv2, latFOM_tray);
#ifdef _OPENMP
				#pragma omp critical(monoP)
#endif
				{
					lattice_result_mono_P.insert(lattice_result_mono_P.end(), latFOM_tray.begin(), latFOM_tray.end());
				}
			}
			if( JudgeSymmetry[Orthorhombic_P] )
			{
				latFOM.putLatticesOfHigherSymmetry(D2h, m_cv2, latFOM_tray);
#ifdef _OPENMP
				#pragma omp critical (ortho_P)
#endif
				{
					lattice_result_ortho_P.insert(lattice_result_ortho_P.end(), latFOM_tray.begin(), latFOM_tray.end());
				}
			}
	   	}
	}

	sort( lattice_result_tri.begin(), lattice_result_tri.end(), VCLatticeFigureOfMeritToCheckSymmetry::cmpFOMdeWolff );
	if( MAX_SIZE < (Int4)lattice_result_tri.size() )
	{
		lattice_result_tri.erase( lattice_result_tri.begin() + MAX_SIZE, lattice_result_tri.end() );
	}

	const size_t num_tri2 = lattice_result_tri.size();
	for(size_t n=0; n<num_tri2; n++)
   	{
		lattice_result_tri[n].setLabel(n+1);
   	}

ZLOG_INFO( "The program has selected "	+ num2str<Int4>( lattice_result_tri.size() )
				+ " triclinic solutions using the Wu figure of merit.\n\n" );


/*   2011.10.19 VIC Tamura */
CHECK_INTERRUPTION();

/*   2011.10.19 VIC Tamura */
LOOP_COUNTER = 0;
Int4 SUM = 0;
for(Int4 i=0; i<NUM_LS; i++) { SUM += lattice_result[i].size(); }

	for(Int4 i=1; i<NUM_LS; i++)
	{
		if( !JudgeSymmetry[i] ) continue;
//		sort( lattice_result[i].begin(), lattice_result[i].end() );

		const Int4 num_lattice = lattice_result[i].size();

#ifdef _OPENMP
		#pragma omp parallel
#endif
		{
			vector< VecDat3<Int4> > closest_hkl_tray;
			vector<bool> is_cal_Q_observed_tray;
			vector<VCLatticeFigureOfMeritToCheckSymmetry> latFOM_tray;

#ifdef _OPENMP
			#pragma omp for
#endif
			for(Int4 index=0; index<num_lattice; index++)
			{
/*   2011.10.19 VIC Tamura */
SET_PROGRESS(99+1*(LOOP_COUNTER++)/SUM, 66, 30); // critical, but works
if(IS_CANSELED()) continue;

				VCLatticeFigureOfMeritToCheckSymmetry& latFOM0 = lattice_result[i][index];
				latFOM0.setLabel(index+1);

				latFOM0.setFigureOfMerit(m_num_ref_figure_of_merit,
											VCData::putPeakQData(),
											closest_hkl_tray, is_cal_Q_observed_tray);

				LatticeFigureOfMeritZeroShift latFOM2 = latFOM0.putLatticeFigureOfMerit();
				pair<bool, ZErrorMessage> ans = latFOM2.fitLatticeParameterLinear(VCData::putPeakQData(),
															closest_hkl_tray, is_cal_Q_observed_tray, false);
				if( ans.first )
				{
					assert( latFOM0.putLatticeFigureOfMerit().putFiguresOfMerit().putNumberOfReflectionsForFigureOfMerit() > 0 );
					latFOM2.setFigureOfMerit(latFOM0.putLatticeFigureOfMerit().putFiguresOfMerit().putNumberOfReflectionsForFigureOfMerit(),
												VCData::putPeakQData(),
												closest_hkl_tray, is_cal_Q_observed_tray);
					if( LatticeFigureOfMerit::cmpFOMWu( latFOM2, latFOM0.putLatticeFigureOfMerit() ) )
					{
						latFOM0.setLatticeFigureOfMerit(latFOM2);
					}
				}

				const LatticeFigureOfMerit::SetOfFigureOfMerit& setFOM = latFOM0.putLatticeFigureOfMerit().putFiguresOfMerit();
				if( setFOM.putFigureOfMeritWu() < MIN_NormM ) continue;
				if( setFOM.putReversedFigureOfMerit() < MIN_RevM ) continue;

				if( eBravaisType(i) == Monoclinic_P )
				{
					if( JudgeSymmetry[Hexagonal] )
					{
						latFOM0.putLatticesOfHigherSymmetry(D6h, m_cv2, latFOM_tray);
#ifdef _OPENMP
						#pragma omp critical (hex)
#endif
						{
							lattice_result_hex.insert(lattice_result_hex.end(), latFOM_tray.begin(), latFOM_tray.end());
						}
					}
				}
				else if( eBravaisType(i) == Monoclinic_B )
				{
					if( JudgeSymmetry[Orthorhombic_C] )
					{
						latFOM0.putLatticesOfHigherSymmetry(D2h, m_cv2, latFOM_tray);
#ifdef _OPENMP
						#pragma omp critical (ortho_B)
#endif
						{
							lattice_result_ortho_B.insert(lattice_result_ortho_B.end(), latFOM_tray.begin(), latFOM_tray.end());
						}
					}
				}
				else if( eBravaisType(i) == Orthorhombic_P )
				{
					if( JudgeSymmetry[Tetragonal_P] )
					{
						latFOM0.putLatticesOfHigherSymmetry(D4h_Z, m_cv2, latFOM_tray);
#ifdef _OPENMP
						#pragma omp critical (tetra_P)
#endif
						{
							lattice_result_tetra_P.insert(lattice_result_tetra_P.end(), latFOM_tray.begin(), latFOM_tray.end());
						}
					}
					if( JudgeSymmetry[Cubic_P] )
					{
						latFOM0.putLatticesOfHigherSymmetry(Oh, m_cv2, latFOM_tray);
#ifdef _OPENMP
						#pragma omp critical (cubic_P)
#endif
						{
							lattice_result_cubic_P.insert(lattice_result_cubic_P.end(), latFOM_tray.begin(), latFOM_tray.end());
						}
					}
				}
				else if( eBravaisType(i) == Orthorhombic_I )
				{
					if( JudgeSymmetry[Tetragonal_I] )
					{
						latFOM0.putLatticesOfHigherSymmetry(D4h_Z, m_cv2, latFOM_tray);
#ifdef _OPENMP
						#pragma omp critical (tetra_I)
#endif
						{
							lattice_result_tetra_I.insert(lattice_result_tetra_I.end(), latFOM_tray.begin(), latFOM_tray.end());
						}
					}
					if( JudgeSymmetry[Cubic_I] )
					{
						latFOM0.putLatticesOfHigherSymmetry(Oh, m_cv2, latFOM_tray);
#ifdef _OPENMP
						#pragma omp critical (cubic_I)
#endif
						{
							lattice_result_cubic_I.insert(lattice_result_cubic_I.end(), latFOM_tray.begin(), latFOM_tray.end());
						}
					}
				}
				else if( eBravaisType(i) == Orthorhombic_F )
				{
					if( JudgeSymmetry[Cubic_F] )
					{
						latFOM0.putLatticesOfHigherSymmetry(Oh, m_cv2, latFOM_tray);
#ifdef _OPENMP
						#pragma omp critical (cubic_F)
#endif
						{
							lattice_result_cubic_F.insert(lattice_result_cubic_F.end(), latFOM_tray.begin(), latFOM_tray.end());
						}
					}
				}
			}
		}
/*   2011.10.19 VIC Tamura */
CHECK_INTERRUPTION();

		sort( lattice_result[i].begin(), lattice_result[i].end(), VCLatticeFigureOfMeritToCheckSymmetry::cmpFOMdeWolff );
		if( MAX_SIZE < (Int4)lattice_result[i].size() )
		{
			lattice_result[i].erase( lattice_result[i].begin() + MAX_SIZE, lattice_result[i].end() );
		}

		const size_t num_lattice2 = lattice_result[i].size();
		for(size_t n=0; n<num_lattice2; n++)
	   	{
			lattice_result[i][n].setLabel(n+1);
	   	}

ZLOG_INFO( "(" + num2str( i+1 ) + ") The number of candidates for " + put_bravais_type_name(eBravaisType(i), abc_axis)
			+ " : " + num2str<Int4>( lattice_result[i].size() ) + "\n" );
	}
ZLOG_INFO( "\n" );
    }
    catch(bad_alloc& ball)
    {
    	throw nerror(ball, __FILE__, __LINE__, __FUNCTION__);
    }
}


void SortingLattice::putLatticeCandidatesForEachBravaisTypes(const vector<SymMat43_VCData>& S_super,
					const Double& MIN_NormM,
					const Double& MIN_RevM,
					const Int4& MAX_SIZE,
					const eABCaxis& abc_axis,
					const eRHaxis& rh_axis,
					vector<VCLatticeFigureOfMeritToCheckSymmetry> lattice_result[NUM_LS]) const
{
	vector<VCLatticeFigureOfMeritToCheckSymmetry>& lattice_result_tri = lattice_result[(size_t)Triclinic];
	putLatticeCandidatesForTriclinic(S_super, MIN_NormM, MIN_RevM, lattice_result_tri);

ZLOG_INFO( "Determination of the Bravais type is being carried out...\n(Solutions with " + putLabel(SCWuM) + " less than " + num2str(MIN_NormM)
		+ " or " + putLabel(SCRevM) + " less than " + num2str(MIN_RevM)
		+ " are automatically removed).\n" );
ZLOG_INFO( "All the unit-cell parameters are being optimized by linear least squares...\n" );

//ZLOG_INFO( "The program has removed "	+ num2str<Int4>( S_super.size() - lattice_result_tri.size() )
//			+ " triclinic solutions because their " + putLabel(SCWuM) + " is less than " + num2str(MIN_NormM)
//			+ " or their " + putLabel(SCRevM) + " is less than " + num2str(MIN_RevM) + ".\n\n" );
//ZLOG_INFO( "Determination of the Bravais type is being carried out...\n" );
	putLatticeCandidatesForEachBravaisTypes(MIN_NormM, MIN_RevM, MAX_SIZE, abc_axis, rh_axis, lattice_result);
}


void SortingLattice::setNumberOfNeighbors(const eABCaxis& baxis_flag,
		bool (*CmpFunc)(const VCLatticeFigureOfMeritToCheckSymmetry&, const VCLatticeFigureOfMeritToCheckSymmetry&),
		vector<VCLatticeFigureOfMeritToCheckSymmetry> lattice_result[NUM_LS]) const
{

#ifdef _OPENMP
		#pragma omp for
#endif
	for(Int4 i=0; i<NUM_LS; i++)
	{
		if( !OutputSymmetry[(size_t)i] ) continue;

		stable_sort( lattice_result[i].begin(), lattice_result[i].end() ); // Sort by the unit-cell volume.
		for(vector<VCLatticeFigureOfMeritToCheckSymmetry>::iterator it=lattice_result[i].begin(); it<lattice_result[i].end(); it++)
		{
			it->putNumberOfLatticesInNeighborhood() = 0;
		}
	}

	const Double coef_lower = 1.0 - m_resol*3.0;
	const Double coef_upper = 1.0 + m_resol*3.0;

	Vec_INT index_tray(put_number_of_bravais_types(), 0);

/*   2011.10.19 VIC Tamura */
Int4 SUM=0, LOOP_COUNTER=0;
for(Int4 i=0; i<NUM_LS; i++ ) { SUM += lattice_result[(size_t)i].size(); }

#ifdef _OPENMP
			#pragma omp for
#endif
	for(Int4 i=0; i<NUM_LS; i++)
	{
		if( !OutputSymmetry[(size_t)i] ) continue;

		const size_t num_lattice = lattice_result[i].size();

		for(size_t index=0; index<num_lattice; index++)
		{
/*   2011.10.19 VIC Tamura */
SET_PROGRESS(100*(LOOP_COUNTER++)/SUM, 97, 1); // critical, but works
if(IS_CANSELED()) continue;

			VCLatticeFigureOfMeritToCheckSymmetry& latFOM0 = lattice_result[i][index];
			const LatticeFigureOfMerit& latFOM0_prim = latFOM0.putLatticeFigureOfMerit();
			if( latFOM0.putNumberOfLatticesInNeighborhood() < 0 ) continue;

			const Double& detS = latFOM0_prim.putDeterminantOfGramMatrix();
			const size_t ibegin = distance( lattice_result[i].begin(), lower_bound( lattice_result[i].begin(), lattice_result[i].end(), detS*coef_lower ) );
			const size_t iend = distance( lattice_result[i].begin(), upper_bound( lattice_result[i].begin()+ibegin, lattice_result[i].end(), detS*coef_upper ) );

			Int4 count=0;
			if( i == (size_t)Triclinic )
			{
				const ReducedLatticeToCheckEquiv RLCS(m_resol, latFOM0_prim.putSellingReducedForm());
				for(size_t index2=ibegin; index2<iend; index2++)
				{
					if( index2 == index ) continue;

					VCLatticeFigureOfMeritToCheckSymmetry& latFOM2 = lattice_result[i][index2];
					const LatticeFigureOfMerit& latFOM2_prim = latFOM2.putLatticeFigureOfMerit();

					// lattice_result_tri[index2] equals trans_mat * RLCB.m_S_super_obtuse * Transpose(trans_mat)
					if( RLCS.equiv( latFOM2_prim.putSellingReducedForm() ) )
					{
						// Compare the figures of merit.
						if( CmpFunc( latFOM2, latFOM0 ) )
						{
							if( latFOM2.putNumberOfLatticesInNeighborhood() >= 0 )
							{
								latFOM2.putNumberOfLatticesInNeighborhood() += 1;
								count = -1;
								break;
							}
						}
						else
						{
							count++;
							latFOM2.putNumberOfLatticesInNeighborhood() = -1;
						}
					}
				}
			}
			else
			{
				for(size_t index2=ibegin; index2<iend; index2++)
				{
					if( index2 == index ) continue;

					VCLatticeFigureOfMeritToCheckSymmetry& latFOM2 = lattice_result[i][index2];
					const LatticeFigureOfMerit& latFOM2_prim = latFOM2.putLatticeFigureOfMerit();

					// *it2 equals trans_mat * RLCS.m_S_super_obtuse * Transpose(trans_mat)
					if( check_equiv_m( latFOM0_prim.putInverseOfBuergerReducedForm(), latFOM2_prim.putInverseOfBuergerReducedForm(), m_resol ) )
					{
						// Compare the figures of merit.
						if( CmpFunc( latFOM2, latFOM0 ) )
						{
							if( latFOM2.putNumberOfLatticesInNeighborhood() >= 0 )
							{
								latFOM2.putNumberOfLatticesInNeighborhood() += 1;
								count = -1;
								break;
							}
						}
						else
						{
							count++;
							latFOM2.putNumberOfLatticesInNeighborhood() = -1;
						}
					}
				}
			}

			latFOM0.putNumberOfLatticesInNeighborhood() = count;
		}

		Int4& index = index_tray[i];
		index = 0;
		for(vector<VCLatticeFigureOfMeritToCheckSymmetry>::const_iterator it=lattice_result[i].begin(); it<lattice_result[i].end(); it++)
	   	{
			if( it->putNumberOfLatticesInNeighborhood() >= 0 ) index++;
	   	}
	}

/*   2011.10.19 VIC Tamura */
CHECK_INTERRUPTION();

for(Int4 i=0; i<NUM_LS; i++)
{
	if( !OutputSymmetry[(size_t)i] ) continue;
	ZLOG_INFO( "(" + num2str( i+1 ) + ") The number of candidates for " + put_bravais_type_name(eBravaisType(i), baxis_flag)
					+ " : " + num2str( lattice_result[i].size() ) + " -> " + num2str( index_tray[i] ) + "\n" );
}
ZLOG_INFO( "\n" );
}

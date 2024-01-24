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
#include "indexing_func_dim3.hh"
#include "lattice_symmetry/gather_q_of_3D_lattice.hh"
#include "utility_data_structure/Bud2.hh"
#include "utility_data_structure/Node3.hh"
#include "utility_data_structure/nrutil_nr.hh"
#include "utility_lattice_reduction/put_Selling_reduced_lattice.hh"
#include "utility_lattice_reduction/put_Buerger_reduced_lattice.hh"
#include "utility_func/lattice_constant.hh"
#include "utility_func/stopx.hh"
#include "utility_func/zmath.hh"
#include "zlog/zlog.hh"
#include "lattice_symmetry/VCLatticeFigureOfMeritToCheckSymmetry.hh"


inline bool operator<(const pair<Double, SymMat43_VCData>& lhs, const pair<Double, SymMat43_VCData>& rhs)
{
	if( lhs.second.first == rhs.second.first ) return false;
	if( lhs.first < rhs.first  ) return true;
	if( lhs.first > rhs.first ) return false;
	return lhs.second.first < rhs.second.first;
}


static bool check_lattice_distance(const SymMat<Double>& S,
		const Double& Min_LatDist)
{
	SymMat<Double> S_super(4);
	NRMat<Int4> trans_mat(4,3);
	if( !put_Selling_reduced_dim_less_than_4(S, S_super, trans_mat, Min_LatDist) ) return false;

	if( S_super(0,0) + S_super(1,1) + S_super(0,1)*2.0 < Min_LatDist ) return false;
	if( S_super(0,0) + S_super(2,2) + S_super(0,2)*2.0 < Min_LatDist ) return false;
	if( S_super(1,1) + S_super(2,2) + S_super(1,2)*2.0 < Min_LatDist ) return false;

	return true;
}



static SymMat43_VCData putTransformMatrixFromSellingReducedToBuergerReduced(
		const SymMat<VCData>& S_super_obtuse)
{
	NRMat<Int4> trans_mat(4,3);

	// S_red = trans_mat * S_super_obtuse * transpose(trans_mat).
	SymMat<VCData> S_red(3);
	putBuergerReducedMatrix(S_super_obtuse, false, S_red, trans_mat);

	return SymMat43_VCData( S_red, put_transform_matrix_row3to4( Inverse3( trans_mat ) ) );
}



// On input, VCData::m_peak_pos[data_type] and bud_in_tree are supposed to be sorted into ascending order.
void combinate_bud(const vector<Bud2>& bud_basket,
		const Double& cv2,
		const Double& min_unit_cell_volume,
		const Double& max_unit_cell_volume,
		const Double& min_length,
	    const Int4& num_ref_figure_of_merit,
	    const Double& MIN_WuM,
	    const ePeakShiftFunctionType& etype_peak_shift,
	    const Double& WlengthX,
	    const vector<ZParawError>& peak_shift_param_rad,
		const Int4& MAX_NODE_NUM2,
		const eConographAnalysisMode& SearchLevel,
		set< pair<Double, SymMat43_VCData> >& node_basket_vol,
		set< pair<Double, SymMat43_VCData> >& node_basket_FOM)
{
	node_basket_vol.clear();
	node_basket_FOM.clear();
	
	const vector<QData>& qdata = VCData::putPeakQData();
	if( (Int4)qdata.size() < 2 || MAX_NODE_NUM2 <= 0 ) return;

	assert(min_unit_cell_volume>=0.0);
	assert(max_unit_cell_volume>0.0);
	const Double Max_CoUNIT_CELL_VOL = 1.0/min_unit_cell_volume;
	const Double Min_CoUNIT_CELL_VOL = 1.0/max_unit_cell_volume;

	const Double Max_detS = Max_CoUNIT_CELL_VOL*Max_CoUNIT_CELL_VOL;
	const Double Min_detS = Min_CoUNIT_CELL_VOL * Min_CoUNIT_CELL_VOL;
	Double Max_detS2 = Max_detS;
	Double MIN_WuM2 = MIN_WuM;

ZLOG_INFO( "Permitted range of primitive-cell volume = " + num2str( 1.0/Max_CoUNIT_CELL_VOL ) + "--" + num2str( 1.0/Min_CoUNIT_CELL_VOL ) + "\n" );
	
	SetSignal(SIGINT);

/*   2011.10.19 VIC Tamura */
Int4 LOOP_COUNTER = 0;

	const Int4 BUDSIZE = bud_basket.size();

#ifdef _OPENMP
		#pragma omp parallel
#endif
	{
		Node3 nodex;
		vector<Bud2>::const_iterator it13_end;

		SymMat<VCData> S(3);
		SymMat<VCData> vc_S_super(4);
		NRMat<Int4> trans_mat(4,3);
		SymMat43_VCData S_red(S, trans_mat);
		bool removed_last_flag;

#ifdef _OPENMP
		#pragma omp for
#endif
		for(Int4 i_=0; i_<BUDSIZE; i_++)
		{
/*   2011.10.19 VIC Tamura */
SET_PROGRESS(100*(LOOP_COUNTER++)/(Int4)bud_basket.size(), 8, 32); // critical, but works
			
			if( putInterruptionSignal() )
			{
				continue;
			}
			const vector<Bud2>::const_iterator it12=bud_basket.begin()+i_;

			const Int4 iK1 = it12->iK1();
			const VCData Q1 = it12->Q1();
			const VCData& Q2 = it12->Q2();
			const VCData& Q12 = it12->Q12();

			it13_end = upper_bound(it12, bud_basket.end(), Bud2::putBud2(iK1+1));

			for(vector<Bud2>::const_iterator it13=it12; it13<it13_end; it13++)
			{
				if( iK1 < 0 )
				{
					if( !vc_equiv(Q1, it13->Q1(), cv2) ) continue;
				}

				const Double InvA = 1.0/Q1.Value();
				const Double s12 = (Q12.Value()-Q1.Value()-Q2.Value())*0.5;

				const VCData& Q3 = it13->Q2();
				const VCData& Q13 = it13->Q12();

				const Double s13 = (Q13.Value()-Q1.Value()-Q3.Value())*0.5;
				const Double B_A = s12*s13*InvA;
				const Double C_A = (Q2.Value()*s13*s13 + Q3.Value()*s12*s12 - Q1.Value()*Q2.Value()*Q3.Value())*InvA;

				// detS = -A*s23^2 + 2*B*s23 - C,
				// A = s11,
				// B = s12*s13,
				// C = s22*s13^2 + s33*s12^2 - s11*s22*s33.
				// (B/A)^2 - C/A - Max_detS/A <= (s23 - B/A)^2 <= (B/A)^2 - C/A - Min_detS/A
				const Double DIFF = sqrt(max(0.0, B_A*B_A - C_A - Min_detS*InvA));
				const Double SumQ = Q12.Value() + Q13.Value() - Q1.Value();
				const Double MinQ123 = SumQ + (B_A - DIFF)*2.0;
				const Double MaxQ123 = SumQ + (B_A + DIFF)*2.0;
				const Int4 Min_iK123 = distance(qdata.begin(), lower_bound(qdata.begin(), qdata.end(), QData(MinQ123, 0.0)));
				const Int4 Max_iK123 = distance(qdata.begin(), upper_bound(qdata.begin()+Min_iK123, qdata.end(), QData(MaxQ123, 0.0)));

				const Double DIFF2 = sqrt(max(0.0, B_A*B_A - C_A - (SearchLevel==ConographQuickSearch?Max_detS2:Max_detS)*InvA));
				const Double MinProhibitedRange = SumQ + (B_A - DIFF2)*2.0;
				const Double MaxProhibitedRange = SumQ + (B_A + DIFF2)*2.0;

				for(Int4 iK123=Min_iK123; iK123<Max_iK123; iK123++)
				{
					if( MinProhibitedRange < qdata[iK123].q && qdata[iK123].q < MaxProhibitedRange ) continue;
					if( !nodex.setBud0(cv2, Q1, Q2, Q3, VCData(iK123, 1), Q12, Q13, S) ) continue;
					if( SearchLevel == ConographQuickSearch )
					{
						if( nodex.putDeterminantS() > Max_detS2 ) continue;

						if( iK1 < 0 )
						{
							const VCData Q_123 = S(0,0) + S(1,1) + S(2,2) + ( S(1,2) - S(0,1) - S(0,2) )*2;
							const pair<vector<QData>::const_iterator, vector<QData>::const_iterator> itpair_qdata = equal_range_qdata(qdata.begin(), qdata.end(), Q_123.Value(), sqrt(Q_123.Variance()*cv2) );
							if( itpair_qdata.first >= itpair_qdata.second ) continue;
						}
					}

					const SymMat<Double> dbl_S = chToDouble(S);
					if( !check_lattice_distance(Inverse3(dbl_S), min_length) ) continue;

					if( !put_Selling_reduced_dim_less_than_4(S, vc_S_super, trans_mat) ) continue;
					moveSmallerDiagonalLeftUpper(vc_S_super);
					S_red = putTransformMatrixFromSellingReducedToBuergerReduced(vc_S_super);

					if( SearchLevel == ConographQuickSearch )
					{
#ifdef _OPENMP
						#pragma omp critical (VOL)
#endif
						{
							node_basket_vol.insert(pair<Double, SymMat43_VCData>(nodex.putDeterminantS(), S_red) );
							if( (Int4)node_basket_vol.size() > MAX_NODE_NUM2 )
							{
								node_basket_vol.erase(--(node_basket_vol.end()));
								Max_detS2 = node_basket_vol.rbegin()->first;
							}
						}
						continue;
					}
					else if( nodex.putDeterminantS() <= Max_detS2 )
					{
						removed_last_flag = false;
#ifdef _OPENMP
						#pragma omp critical (VOL)
#endif
						{
							node_basket_vol.insert(pair<Double, SymMat43_VCData>(nodex.putDeterminantS(), S_red) );
							if( (Int4)node_basket_vol.size() > MAX_NODE_NUM2 )
							{
								S_red = node_basket_vol.rbegin()->second;
								removed_last_flag = true;
								node_basket_vol.erase(--(node_basket_vol.end()));
								Max_detS2 = node_basket_vol.rbegin()->first;
							}
						}
						if( !removed_last_flag ) continue;
					}

					VCLatticeFigureOfMeritToCheckSymmetry latFOM(BravaisType( pair<eCentringType, ePointGroup>(Prim, Ci) ), S_red,
																etype_peak_shift, WlengthX, peak_shift_param_rad);

					// Calculate figures of merit as triclinic
					latFOM.setWuFigureOfMerit( num_ref_figure_of_merit, VCData::putPeakQData(), 0.8, 12.5 );

					const LatticeFigureOfMerit::SetOfFigureOfMerit& setFOM = latFOM.putLatticeFigureOfMerit().putFiguresOfMerit();
					if( setFOM.putFigureOfMeritWu() < MIN_WuM2 ) continue;

#ifdef _OPENMP
					#pragma omp critical (FOM)
#endif
					{
						node_basket_FOM.insert(pair<Double, SymMat43_VCData>( -setFOM.putFigureOfMeritWu(), S_red ) );
						if( (Int4)node_basket_FOM.size() > MAX_NODE_NUM2 )
						{
							node_basket_FOM.erase(--(node_basket_FOM.end()));
							MIN_WuM2 = -node_basket_FOM.rbegin()->first;
						}
					}
				}
			}
		}
	}

/*   2011.10.19 VIC Tamura */
CHECK_INTERRUPTION();
}

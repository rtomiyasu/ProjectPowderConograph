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
#include "../utility_data_structure/FracMat.hh"
#include "../utility_func/chToDouble.hh"
#include "../utility_lattice_reduction/put_Selling_reduced_lattice.hh"
#include "../utility_lattice_reduction/put_Buerger_reduced_lattice.hh"
#include "../laue_group/LaueGroup.hh"
#include "../point_group/PGNormalSeriesTray.hh"
#include "../model_function/LatticeDistanceModel.hh"
#include "../zlog/zlog.hh"
#include "gather_q_of_3D_lattice.hh"
#include "gather_q_of_Ndim_lattice.hh"
#include "ReducedLatticeToCheckBravais.hh"
#include "LatticeFigureOfMerit.hh"
#include "../qc/gather_qcal2.hh"

const Double LatticeFigureOfMerit::m_cv2 = 9.0;

const NRMat<Int4> LatticeFigureOfMerit::m_tmat_prim_to_face = put_transform_matrix_row3to4( transpose( BravaisType::putTransformMatrixFromPrimitiveToFace() ) );
const NRMat<Int4> LatticeFigureOfMerit::m_tmat_prim_to_body = put_transform_matrix_row3to4( BravaisType::putTransformMatrixFromBodyToPrimitive() );
const NRMat<Int4> LatticeFigureOfMerit::m_tmat_prim_to_rhomhex = put_transform_matrix_row3to4( transpose( BravaisType::putTransformMatrixFromPrimitiveToRhomHex() ) );
const NRMat<Int4> LatticeFigureOfMerit::m_tmat_prim_to_base[3] =
		{
				put_transform_matrix_row3to4( transpose( BravaisType::putTransformMatrixFromPrimitiveToBase(BaseA_Axis) ) ),
				put_transform_matrix_row3to4( transpose( BravaisType::putTransformMatrixFromPrimitiveToBase(BaseB_Axis) ) ),
				put_transform_matrix_row3to4( transpose( BravaisType::putTransformMatrixFromPrimitiveToBase(BaseC_Axis) ) )
		};
const NRMat<Int4> LatticeFigureOfMerit::m_tmat_prim_to_prim = put_transform_matrix_row3to4();

LatticeFigureOfMerit::LatticeFigureOfMerit()
	: m_S_optimized( SymMat43_Double( SymMat<Double>(3), NRMat<Int4>(4,3) ) ), m_S_red(3),
	  m_determ_S_red(0.0)
{
}


LatticeFigureOfMerit::LatticeFigureOfMerit(const Double& rhs)
	: m_S_optimized( SymMat43_Double( SymMat<Double>(3), NRMat<Int4>(4,3) ) ), m_S_red(3),
	  m_determ_S_red(rhs)
{
}


LatticeFigureOfMerit::LatticeFigureOfMerit(const BravaisType& brat,
		const SymMat43_Double& S)
	: m_S_optimized( SymMat43_Double( SymMat<Double>(3), NRMat<Int4>(4,3) ) ), m_S_red(3)
{
	this->setLatticeConstants43(brat, S);
}

#ifdef DEBUG
static bool checkInitialLatticeParameters(
		const BravaisType& brat,
		const SymMat<Double>& S_red)
{
	const SymMat<Double> inv_S_red( Inverse3(S_red) );

	if( brat.enumLaueGroup() == C2h_Y && brat.enumCentringType() == Prim )
	{
		assert( inv_S_red(0,2) <= 0.0 &&
				inv_S_red(0,0) * 0.9999 < inv_S_red(2,2)
				&& fabs( inv_S_red(0,2) ) * 1.9999 < inv_S_red(2,2)
				&& fabs( inv_S_red(0,2) ) * 1.9999 < inv_S_red(0,0) );
	}
	else if( brat.enumLaueGroup() == C2h_Z && brat.enumCentringType() == Prim )
	{
		assert( inv_S_red(0,1) <= 0.0
				&& inv_S_red(0,0) * 0.9999 < inv_S_red(1,1)
				&& fabs( inv_S_red(0,1) ) * 1.9999 < inv_S_red(0,0)
				&& fabs( inv_S_red(0,1) ) * 1.9999 < inv_S_red(1,1) );
	}
	else if( brat.enumLaueGroup() == C2h_X && brat.enumCentringType() == Prim )
	{
		assert( inv_S_red(1,2) <= 0.0
				&& inv_S_red(1,1) * 0.9999 < inv_S_red(2,2)
				&& fabs( inv_S_red(1,2) ) * 1.9999 < inv_S_red(1,1)
				&& fabs( inv_S_red(1,2) ) * 1.9999 < inv_S_red(2,2) );
	}
	else if( brat.enumLaueGroup() == C2h_Y && brat.enumCentringType() == BaseZ )
	{
		assert( inv_S_red(0,2) <= 0.0
				&& fabs( inv_S_red(0,2) ) * 0.9999 < inv_S_red(2,2)
				&& fabs( inv_S_red(0,2) ) * 1.9999 < inv_S_red(0,0) );
	}
	else if( brat.enumLaueGroup() == C2h_Z && brat.enumCentringType() == BaseX )
	{
		assert( inv_S_red(0,1) <= 0.0
				&& fabs( inv_S_red(0,1) ) * 0.9999 < inv_S_red(0,0)
				&& fabs( inv_S_red(0,1) ) * 1.9999 < inv_S_red(1,1) );
	}
	else if( brat.enumLaueGroup() == C2h_X && brat.enumCentringType() == BaseY )
	{
		assert( inv_S_red(1,2) <= 0.0
				&& fabs( inv_S_red(1,2) ) * 0.9999 < inv_S_red(1,1)
				&& fabs( inv_S_red(1,2) ) * 1.9999 < inv_S_red(2,2) );
	}
	else if( brat.enumBravaisType() == Orthorhombic_C )
	{
		assert( brat.enumCentringType() == BaseZ );
		assert( inv_S_red(0,0) * 0.9999 < inv_S_red(1,1) );
	}
	else if( brat.enumLaueGroup() == D2h && brat.enumCentringType() == Prim )
	{
		assert( inv_S_red(0,0) * 0.9999 < inv_S_red(1,1)
				&& inv_S_red(1,1) * 0.9999 < inv_S_red(2,2) );
	}
	return true;
}
#endif

void putTransformMatrixToBuergerReduced(const SymMat<Double>& S, NRMat<Int4>& trans_mat)
{
	assert( S.size() == 3 );

	SymMat<Double> S_super_obtuse(4);
	put_Selling_reduced_dim_less_than_4(S, S_super_obtuse, trans_mat);
	moveSmallerDiagonalLeftUpper(S_super_obtuse, trans_mat);

	// S_red = trans_mat * S_super_obtuse * transpose(trans_mat).
	SymMat<Double> S_red(3);
	NRMat<Int4> trans_mat2;
	putBuergerReducedMatrix(S_super_obtuse, S_red, trans_mat2);
	trans_mat = mprod( trans_mat2, put_transform_matrix_row4to3(trans_mat) );
}


void LatticeFigureOfMerit::setInverseOfBuergerReducedForm(NRMat<Int4>& trans_mat)
{
	if( m_brat.enumBravaisType() == Triclinic )
	{
		// trans_mat * Inverse(m_S_optimized.first) * transpose(trans_mat) is Buerger-reduced
		// <=> Inverse of transpose(Inverse(trans_mat)) * m_S_optimized.first * Inverse(trans_mat) is Buerger-reduced.
		putTransformMatrixToBuergerReduced(Inverse3(m_S_optimized.first), trans_mat);
		transpose_square_matrix(trans_mat);
		m_S_red = transform_sym_matrix(Inverse3(trans_mat), m_S_optimized.first);
	}
	else
	{
		m_S_red = m_S_optimized.first;
		trans_mat = identity_matrix<Int4>(3);
		if( m_brat.enumBravaisType() == Monoclinic_P )
		{
			if( m_brat.enumLaueGroup() == C2h_X )
			{
				putBuergerReducedMonoclinicP(1, 2, m_S_red, trans_mat);
			}
			else if( m_brat.enumLaueGroup() == C2h_Y )
			{
				putBuergerReducedMonoclinicP(0, 2, m_S_red, trans_mat);
			}
			else //if( m_brat.enumLaueGroup() == C2h_Z )
			{
				putBuergerReducedMonoclinicP(0, 1, m_S_red, trans_mat);
			}
		}
		else if( m_brat.enumBravaisType() == Monoclinic_B )
		{
			m_S_red = m_S_optimized.first;
			putBuergerReducedMonoclinicB(m_brat, m_S_red, trans_mat);
		}
		else if( m_brat.enumLaueGroup() == D2h )
		{
			m_S_red = m_S_optimized.first;
			putBuergerReducedOrthorhombic(m_brat.enumCentringType(), m_S_red, trans_mat);
		}
	}

	assert( checkInitialLatticeParameters(m_brat, m_S_red) );
}


// This method assumes that S.second * S.first * Transpose(S.second) is obtuse.
void LatticeFigureOfMerit::setLatticeConstants43(const BravaisType& brat, const SymMat43_Double& S)
{
	m_brat = brat;
	m_S_optimized = S;

	NRMat<Int4> trans_mat;
	setInverseOfBuergerReducedForm(trans_mat);	// Set m_S_red from m_S_optimized.

	m_determ_S_red = Determinant3( m_S_optimized.first );
	m_figures_of_merit.reset();
}


ZErrorMessage LatticeFigureOfMerit::setLatticeConstants(const BravaisType& brat, const SymMat<Double>& Sval)
{
	assert( Sval.size()==3 );

	SymMat43_Double S_red_optimized = SymMat43_Double(Sval, NRMat<Int4>(4,3));
	cal_average_crystal_system(brat.enumLaueGroup(), S_red_optimized.first);
	if( brat.enumCentringType() == Face )
	{
		S_red_optimized.second = m_tmat_prim_to_face;
	}
	else if( brat.enumCentringType() == Inner )
	{
		S_red_optimized.second = m_tmat_prim_to_body;
	}
	else if( brat.enumCentringType() == BaseX
			|| brat.enumCentringType() == BaseY
			|| brat.enumCentringType() == BaseZ )
	{
		S_red_optimized.second = m_tmat_prim_to_base[ (size_t)brat.enumBASEaxis() ];
	}
	else if( brat.enumCentringType() == Rhom_hex )
	{
		S_red_optimized.second = m_tmat_prim_to_rhomhex;
	}
	else // if( brat.enumCentringType() == Prim )
	{
		S_red_optimized.second = m_tmat_prim_to_prim;
	}

	// S_super_obtuse = trans_mat * S_red.first * Transpose(trans_mat).
	SymMat<Double> S_super_obtuse = transform_sym_matrix(S_red_optimized.second, S_red_optimized.first);
	if( !put_Selling_reduced_dim_less_than_4(S_super_obtuse, S_red_optimized.second) )
	{
		return ZErrorMessage(ZErrorArgument, "The argument matrix is not positive definite" __FILE__, __LINE__, __FUNCTION__);
	}
	moveSmallerDiagonalLeftUpper(S_super_obtuse, S_red_optimized.second);

	setLatticeConstants43(brat, S_red_optimized);

	return ZErrorMessage();
}


inline bool checkIfFirstEntryIsPositive(const VecDat3<Int4>& rhs)
{
	for(Int4 i=0; i<3; i++)
	{
		if( rhs[i] == 0 ) continue;
		if( rhs[i] > 0 ) return true;
		else return false;
	}
	return false;
}


static bool cmp_func(const HKL_Q& lhs, const HKL_Q& rhs)
{
	if( lhs.Q() < rhs.Q() ) return true;
	else if( lhs.Q() > rhs.Q() ) return false;
	return abs(lhs.HKL()[0])+abs(lhs.HKL()[1])+abs(lhs.HKL()[2]) < abs(rhs.HKL()[0])+abs(rhs.HKL()[1])+abs(rhs.HKL()[2]);
};


void LatticeFigureOfMerit::putMillerIndicesInRange(const Double& qend,
		const Int4& irc_type,
		vector<HKL_Q>& cal_hkl_tray) const
{
	cal_hkl_tray.clear();

	vector<HKL_Q> cal_hkl_tray2;
	gatherQcal(this->putSellingReducedForm(), qend, m_S_optimized.second, this->putBravaisType(), irc_type, cal_hkl_tray2);

	set< VecDat3<Int4> > found_hkl_tray;
	vector<MillerIndex3> equiv_hkl_tray;
	
	PGNormalSeriesTray normal_tray(m_brat.enumLaueGroup());
	LaueGroup lg(m_brat.enumLaueGroup());
	VecDat3<Int4> hkl;
	
	for(vector<HKL_Q>::const_iterator it=upper_bound(cal_hkl_tray2.begin(), cal_hkl_tray2.end(), HKL_Q(NRVec<Int4>(), 0.0));
			it<cal_hkl_tray2.end(); it++)
	{
		hkl = it->HKL();
		if( found_hkl_tray.find(hkl) != found_hkl_tray.end() )
		{
			continue;
		}
		if( !checkIfFirstEntryIsPositive(hkl) ) hkl *= -1;
		
		normal_tray.putHKLEquiv(MillerIndex3(hkl[0], hkl[1], hkl[2]), equiv_hkl_tray);
#ifdef DEBUG
		if( (Int4)equiv_hkl_tray.size() != lg->LaueMultiplicity(hkl) )
		{
ZLOG_INFO( num2str(hkl[0]) + " "
		+ num2str(hkl[1]) + " "
		+ num2str(hkl[2]) + "\n"
		+ num2str( equiv_hkl_tray.size() ) + "\n"
		+ num2str( lg->LaueMultiplicity(hkl) ) + "\n" );
		}
#endif

		for(vector<MillerIndex3>::const_iterator ithkl=equiv_hkl_tray.begin(); ithkl<equiv_hkl_tray.end(); ithkl++)
		{
			found_hkl_tray.insert( VecDat3<Int4>( (*ithkl)[0], (*ithkl)[1], (*ithkl)[2] ) );
		}
		
		cal_hkl_tray.push_back( HKL_Q(hkl, it->Q()) );
	}
	sort( cal_hkl_tray.begin(), cal_hkl_tray.end(), cmp_func );
}


void LatticeFigureOfMerit::setFigureOfMerit(const Int4& num_ref_figure_of_merit,
		const vector<QData>& qdata,
		vector< VecDat3<Int4> >& closest_hkl_tray,
		Vec_BOOL& Q_observed_flag)
{
	assert( num_ref_figure_of_merit <= (Int4)qdata.size() );

	// Qdata is sorted into ascended order.
	m_figures_of_merit.reset();
	m_figures_of_merit.putNumberOfReflectionsForFigureOfMerit() = num_ref_figure_of_merit;

	const Int4& num_Q = m_figures_of_merit.putNumberOfReflectionsForFigureOfMerit();
	closest_hkl_tray.clear();
	Q_observed_flag.clear();
	closest_hkl_tray.resize(num_Q, 0);
	Q_observed_flag.resize(num_Q, false);

	if( num_Q <= 0 ) return;

//	const Double MinQ = qdata[0].q - sqrt( m_cv2 * qdata[0].q_var );
	const Double MaxQ = qdata[num_Q-1].q + sqrt( m_cv2 * qdata[num_Q-1].q_var );
	const SymMat<Double> S_sup( this->putSellingReducedForm() );

	vector<HKL_Q> cal_hkl_tray;
	gatherQcal(S_sup, MaxQ, m_S_optimized.second, cal_hkl_tray);
	if( cal_hkl_tray.empty() ) return;

	vector< vector<HKL_Q>::const_iterator > closest_hkl_q_tray;
	associateQobsWithQcal(qdata.begin(), qdata.begin()+num_Q, cal_hkl_tray, closest_hkl_q_tray);
	const vector<HKL_Q>::const_iterator it_begin = closest_hkl_q_tray[0];
	const vector<HKL_Q>::const_iterator it_end = closest_hkl_q_tray[num_Q-1] + 1;
	assert( it_end <= cal_hkl_tray.end() );
	if( it_begin + 1 >= it_end ) return;

	Double diff;
	Double actually_disc = 0.0;
	Int4 num_q_observed = 0;
	for(Int4 k=0; k<num_Q; k++)
	{
		closest_hkl_tray[k] = closest_hkl_q_tray[k]->HKL();
		diff = qdata[k].q - closest_hkl_q_tray[k]->Q();
		actually_disc += fabs( diff );
		if( diff * diff <= m_cv2 * qdata[k].q_var )
		{
			Q_observed_flag[k] = true;
			num_q_observed++;
		}
		else Q_observed_flag[k] = false;
	}
	actually_disc /= num_Q;
	m_figures_of_merit.putNumQobsAssociatedWithCloseHKL() = num_q_observed;

	vector< vector<QData>::const_iterator > closest_q_tray;
	associateQcalWithQobs(it_begin, it_end, qdata, closest_q_tray);

	const LaueGroup lg(m_brat.enumLaueGroup());

	Double inv_mult = 2.0 / lg->LaueMultiplicity( it_begin->HKL() );
	Double num_total_hkl = inv_mult;
	Double rev_actually_disc = fabs( it_begin->Q() - closest_q_tray[0]->q ) * inv_mult;

	Double sum_diff = 0.0;
	Int4 index = 1;
	for(vector<HKL_Q>::const_iterator it=it_begin+1; it<it_end; it++, index++)
	{
		inv_mult = 2.0 / lg->LaueMultiplicity( it->HKL() );
		num_total_hkl += inv_mult;
		rev_actually_disc += fabs( it->Q() - closest_q_tray[index]->q ) * inv_mult;

		diff = it->Q() - (it-1)->Q();
		sum_diff += diff * diff;
	}
	m_figures_of_merit.putContinuousNumberOfHKLInRange() = num_total_hkl;
	rev_actually_disc /= num_total_hkl;

	// Calculate the symmetric figures of merit by Wolff.
	m_figures_of_merit.putFigureOfMeritWolff() = ( (it_end - 1)->Q() - it_begin->Q() ) / ( 2.0*actually_disc*num_total_hkl );
	m_figures_of_merit.putFigureOfMeritWu() = sum_diff / ( 4.0 * actually_disc * ( (it_end - 1)->Q() - it_begin->Q() ) );
	m_figures_of_merit.putReversedFigureOfMerit() = ( qdata[num_Q-1].q - qdata[0].q ) / ( 2.0*rev_actually_disc*num_Q );
}




void LatticeFigureOfMerit::setDeWolffFigureOfMerit(
		const Int4& num_ref_figure_of_merit,
		const Int4& irc_type,
		const vector<QData>& qdata)
{
	assert( num_ref_figure_of_merit <= (Int4)qdata.size() );

	// Qdata is sorted into ascended order.
	m_figures_of_merit.reset();
	m_figures_of_merit.putNumberOfReflectionsForFigureOfMerit() = num_ref_figure_of_merit;

	const Int4& num_Q = m_figures_of_merit.putNumberOfReflectionsForFigureOfMerit();
	if( num_Q <= 0 ) return;

	const Double MaxQ = qdata[num_Q-1].q + sqrt( m_cv2 * qdata[num_Q-1].q_var );
	const SymMat<Double> S_sup( this->putSellingReducedForm() );

	vector<HKL_Q> cal_hkl_tray;
	gatherQcal(S_sup, MaxQ, m_S_optimized.second, this->putBravaisType(), irc_type, cal_hkl_tray);
	if( cal_hkl_tray.empty() ) return;

	vector< vector<HKL_Q>::const_iterator > closest_hkl_q_tray;
	associateQobsWithQcal(qdata.begin(), qdata.begin()+num_Q, cal_hkl_tray, closest_hkl_q_tray);
	const vector<HKL_Q>::const_iterator it_begin = closest_hkl_q_tray[0];
	const vector<HKL_Q>::const_iterator it_end = closest_hkl_q_tray[num_Q-1] + 1;
	assert( it_end <= cal_hkl_tray.end() );
	if( it_begin + 1 >= it_end ) return;

	Double actually_disc = 0.0;
	for(Int4 k=0; k<num_Q; k++)
	{
		actually_disc += fabs( qdata[k].q - closest_hkl_q_tray[k]->Q() );
	}
	actually_disc /= num_Q;

	const LaueGroup lg(m_brat.enumLaueGroup());
	Double num_total_hkl = 0.0;
	for(vector<HKL_Q>::const_iterator it=it_begin; it<it_end; it++)
	{
		num_total_hkl += 2.0 / lg->LaueMultiplicity( it->HKL() );
	}

	// Calculate the symmetric figures of merit by Wolff.
	m_figures_of_merit.putFigureOfMeritWolff() = ( (it_end - 1)->Q() - it_begin->Q() ) / ( 2.0*actually_disc*num_total_hkl );
}



void LatticeFigureOfMerit::setWuFigureOfMerit(const Int4& num_ref_figure_of_merit,
		const vector<QData>& qdata,
		const Double& min_thred_num_hkl,
		const Double& max_thred_num_hkl)
{
	m_figures_of_merit.reset();
	m_figures_of_merit.putNumberOfReflectionsForFigureOfMerit() = min( num_ref_figure_of_merit, (Int4)qdata.size() );
	const Int4& num_Q = m_figures_of_merit.putNumberOfReflectionsForFigureOfMerit();
	if( num_Q <= 0 ) return;

	const Double MinQ = qdata[0].q - sqrt( m_cv2 * qdata[0].q_var );
	const Double MaxQ = qdata[num_Q-1].q + sqrt( m_cv2 * qdata[num_Q-1].q_var );

	const SymMat<Double> S_sup( this->putSellingReducedForm() );

	vector<HKL_Q> cal_hkl_tray;
	gatherQcal(S_sup, MaxQ, cal_hkl_tray);
	const Double num_hkl_in_range = distance( lower_bound(cal_hkl_tray.begin(), cal_hkl_tray.end(), HKL_Q(0,MinQ)), cal_hkl_tray.end() );
	if( num_hkl_in_range < num_Q * min_thred_num_hkl ) return;
	if( num_hkl_in_range > num_Q * max_thred_num_hkl ) return;

	vector< vector<HKL_Q>::const_iterator > closest_hkl_q_tray;
	associateQobsWithQcal(qdata.begin(), qdata.begin()+num_Q, cal_hkl_tray, closest_hkl_q_tray);
	const vector<HKL_Q>::const_iterator it_begin = closest_hkl_q_tray[0];
	const vector<HKL_Q>::const_iterator it_end = closest_hkl_q_tray[num_Q-1] + 1;
	assert( it_end <= cal_hkl_tray.end() );
	if( it_begin + 1 >= it_end ) return;

	Double actually_disc = 0.0;
	for(Int4 k=0; k<num_Q; k++)
	{
		actually_disc += fabs( qdata[k].q - closest_hkl_q_tray[k]->Q() );
	}
	actually_disc /= num_Q;

	Double sum_diff = 0.0, diff;
	for(vector<HKL_Q>::const_iterator it=it_begin+1; it<it_end; it++)
	{
		diff = it->Q() - (it-1)->Q();
		sum_diff += diff * diff;
	}

	// Calculate the figure of merit by Wu.
	m_figures_of_merit.putFigureOfMeritWu() = sum_diff / ( 4.0 * actually_disc * ( (it_end - 1)->Q() - it_begin->Q() ) );
}


pair<bool, ZErrorMessage> LatticeFigureOfMerit::fitLatticeParameterLinear(
		const vector<QData>& qdata,
		const vector< VecDat3<Int4> >& hkl_to_fit,
		const vector<bool>& fix_or_fit_flag, const bool& output_view_flag)
{
	const size_t isize = hkl_to_fit.size();
	
	assert( hkl_to_fit.size() == fix_or_fit_flag.size() );
	assert( hkl_to_fit.size() <= qdata.size() );

	Vec_DP ydata(isize), ydata_err(isize);
	Vec_BOOL nxfit(isize);
	Int4 data_num=0;
	
	for(size_t i=0; i<isize; i++)
	{
		ydata[i] = qdata[i].q;
		ydata_err[i] = sqrt_d( qdata[i].q_var );
		if( ydata_err[i] <= 0.0 )
		{
			nxfit[i] = false;
		}
		else
		{
			nxfit[i] = fix_or_fit_flag[i];
			if( nxfit[i] ) data_num++;
		}
	}
	
	LaueGroup lg(m_brat.enumLaueGroup());
	Mat_DP_constr cmat;
	lg->putLatticeConstantFlag(cmat);
	if( data_num <= countNumberOfIndependentParam(cmat.begin(),cmat.end()) )
	{
		return pair<bool, ZErrorMessage>(false, ZErrorMessage("NUMBER OF DATA IS TOO SMALL", __FILE__, __LINE__, __FUNCTION__));
	}
	setIndex(cmat);

	vector<Double> init_param(6);
	const SymMat<Double>& S_val = this->putOptimizedForm().first;
	init_param[0] = S_val(0,0);
	init_param[1] = S_val(1,1);
	init_param[2] = S_val(2,2);
	init_param[3] = S_val(1,2);
	init_param[4] = S_val(0,2);
	init_param[5] = S_val(0,1);
	
	LatticeDistanceModel latModel;
	latModel.setConstraint(&cmat[0]);
	Double chisq_all;
	pair<bool, ZErrorMessage> ans = latModel.setFittedParam(hkl_to_fit, ydata, ydata_err, nxfit,
														output_view_flag, 0.0, 1, init_param, chisq_all);
	if( !(ans.first)) return ans;

	LatticeFigureOfMerit new_lat(*this);
	SymMat<Double> S_red_optimized(3);
	latModel.putResult(S_red_optimized);
	new_lat.setLatticeConstants(m_brat, S_red_optimized);
	new_lat.setFigureOfMerit( m_figures_of_merit.putNumberOfReflectionsForFigureOfMerit(), qdata );

	if( cmpFOMdeWolff(new_lat, *this) )
	{
		*this = new_lat;
		return pair<bool, ZErrorMessage>(true, ZErrorMessage());
	}
	else return pair<bool, ZErrorMessage>(false, ZErrorMessage());
}


void LatticeFigureOfMerit::putEquivalentLatticeConstantsDegreeWithOtherCentring(
		const eABCaxis& abc_axis, const eRHaxis& rh_axis, const Double& resol,
		vector< pair< eBravaisType, SymMat<Double> > >& ans) const
{
	ans.clear();

	// Calculate figures of merit as triclinic
	const ReducedLatticeToCheckBravais RLCB(abc_axis, rh_axis, false, resol, this->putSellingReducedForm());
	const SymMat<Double> S_obtuse = this->putSellingReducedForm();

	if( this->enumBravaisType() != Rhombohedral )
	{
		const map< SymMat<Double>, NRMat<Int4> >& S_red_tray = RLCB.checkCentringType(BravaisType(Rhombohedral, abc_axis, rh_axis));
		for(map< SymMat<Double>, NRMat<Int4> >::const_iterator it=S_red_tray.begin(); it!=S_red_tray.end(); it++)
		{
			ans.push_back( pair< eBravaisType, SymMat<Double> >(Rhombohedral, it->first) );
		}
	}
	if( this->enumCentringType() != Face )
	{
		const map< SymMat<Double>, NRMat<Int4> >& S_red_tray = RLCB.checkCentringType(BravaisType(Orthorhombic_F, abc_axis, rh_axis));
		for(map< SymMat<Double>, NRMat<Int4> >::const_iterator it=S_red_tray.begin(); it!=S_red_tray.end(); it++)
		{
			ans.push_back( pair< eBravaisType, SymMat<Double> >(Orthorhombic_F, it->first) );
		}
	}
	if( this->enumCentringType() != Inner )
	{
		const map< SymMat<Double>, NRMat<Int4> >& S_red_tray = RLCB.checkCentringType(BravaisType(Orthorhombic_I, abc_axis, rh_axis));
		for(map< SymMat<Double>, NRMat<Int4> >::const_iterator it=S_red_tray.begin(); it!=S_red_tray.end(); it++)
		{
			ans.push_back( pair< eBravaisType, SymMat<Double> >(Orthorhombic_I, it->first) );
		}
	}
	if( this->enumCentringType() != BaseX && this->enumCentringType() != BaseY && this->enumCentringType() != BaseZ )
	{
		const map< SymMat<Double>, NRMat<Int4> >& S_red_tray = RLCB.checkCentringType(BravaisType(Monoclinic_B, abc_axis, rh_axis));
		for(map< SymMat<Double>, NRMat<Int4> >::const_iterator it=S_red_tray.begin(); it!=S_red_tray.end(); it++)
		{
			ans.push_back( pair< eBravaisType, SymMat<Double> >(Monoclinic_B, it->first) );
		}
	}

	NRMat<Int4> trans_mat;
	SymMat<Double> S_red(3);
	if( this->enumBravaisType() == Rhombohedral || this->enumCentringType() != Prim )
	{
		const SymMat<Double> S_super = put_sym_matrix_sizeNplus1toN(this->putSellingReducedForm());
		putTransformMatrixToBuergerReduced(Inverse3(S_super), trans_mat);
		transpose_square_matrix(trans_mat);
		ans.push_back( pair< eBravaisType, SymMat<Double> >(Triclinic, transform_sym_matrix(Inverse3(trans_mat), S_super) ) );
	}
}


void LatticeFigureOfMerit::putEquivalentLatticeConstantsDegreeWithOtherCentring(
		const eABCaxis& abc_axis, const eRHaxis& rh_axis, const Double& resol,
		vector< pair< eBravaisType, pair< VecDat3<Double>, VecDat3<Double> > > >& ans) const
{
	vector< pair< eBravaisType, SymMat<Double> > > ans0;
	putEquivalentLatticeConstantsDegreeWithOtherCentring(abc_axis, rh_axis, resol, ans0);

	ans.clear();
	ans.resize( ans0.size() );
	vector< pair< eBravaisType, pair< VecDat3<Double>, VecDat3<Double> > > >::iterator it2 = ans.begin();
	for(vector< pair< eBravaisType, SymMat<Double> > >::const_iterator it=ans0.begin(); it<ans0.end(); it++, it2++)
	{
		it2->first = it->first;
		LatticeFigureOfMerit::putLatticeConstantsDegree( BravaisType(it->first, abc_axis, rh_axis), it->second, abc_axis, rh_axis, it2->second.first, it2->second.second );
	}
}


void LatticeFigureOfMerit::printLatticeInformation(
					const eABCaxis& abc_axis,
					const eRHaxis& rh_axis,
					const Double& resol,
					const Int4& label_start0,
					ostream* os) const
{
	Int4 label_start = label_start0;
	os->width(label_start);
  	*os << "" << "<CrystalSystem>";
	os->width(17);
	*os << put_bravais_type_name(this->enumBravaisType(), abc_axis);
  	*os << " </CrystalSystem>\n\n";

	os->width(label_start); *os << "";
  	*os << "<!-- a, b, c(angstrom), alpha, beta, gamma(deg.)-->\n";

	VecDat3<Double> length_axis, angle_axis;
	if( this->enumBravaisType() == Rhombohedral )
	{
		this->putReducedLatticeConstantsDegree(abc_axis, Rho_Axis, length_axis, angle_axis);

		os->width(label_start); *os << "";
		*os << "<ReducedLatticeParameters axis=\"Rhombohedral\">";
	  	os->width(14);
	  	*os << length_axis[0];
	  	os->width(14);
	   	*os << length_axis[1];
	  	os->width(14);
	   	*os << length_axis[2];
	 	os->width(14);
	   	*os << angle_axis[0];
	  	os->width(14);
	   	*os << angle_axis[1];
	  	os->width(14);
	   	*os << angle_axis [2];
	  	*os << " </ReducedLatticeParameters>\n";

		this->putReducedLatticeConstantsDegree(abc_axis, Hex_Axis, length_axis, angle_axis);

		os->width(label_start); *os << "";
	  	*os << "<ReducedLatticeParameters axis=\"Hexagonal\">";
	  	os->width(14);
	  	*os << length_axis[0];
	  	os->width(14);
	   	*os << length_axis[1];
	  	os->width(14);
	   	*os << length_axis[2];
	 	os->width(14);
	   	*os << angle_axis[0];
	  	os->width(14);
	   	*os << angle_axis[1];
	  	os->width(14);
	   	*os << angle_axis[2];
	  	*os << " </ReducedLatticeParameters>\n\n";
	}
	else
	{
		this->putReducedLatticeConstantsDegree(abc_axis, Rho_Axis, length_axis, angle_axis);

		os->width(label_start); *os << "";
		*os << "<ReducedLatticeParameters>";
	  	os->width(14);
	  	*os << length_axis[0];
	  	os->width(14);
	   	*os << length_axis[1];
	  	os->width(14);
	   	*os << length_axis[2];
	 	os->width(14);
	   	*os << angle_axis[0];
	  	os->width(14);
	   	*os << angle_axis[1];
	  	os->width(14);
	   	*os << angle_axis[2];
	  	*os << " </ReducedLatticeParameters>\n";
	}

	this->putOptimizedLatticeConstantsDegree(abc_axis, rh_axis, length_axis, angle_axis);

	os->width(label_start); *os << "";
	*os << "<OptimizedLatticeParameters>";
  	os->width(14);
  	*os << length_axis[0];
  	os->width(14);
   	*os << length_axis[1];
  	os->width(14);
   	*os << length_axis[2];
 	os->width(14);
   	*os << angle_axis[0];
  	os->width(14);
   	*os << angle_axis[1];
  	os->width(14);
   	*os << angle_axis[2];
  	*os << " </OptimizedLatticeParameters>\n\n";

	os->width(label_start); *os << "";
  	if( this->enumBravaisType() == Rhombohedral )
  	{
		if( rh_axis == Hex_Axis )
		{
			*os << "<VolumeOfUnitCell axis=\"Hexagonal\">";
	  	  	os->width(14);
	  	  	*os << this->putLatticeVolume();
		}
		else // if( rh_axis == Rho_Axis )
		{
			*os << "<VolumeOfUnitCell axis=\"Rhombohedral\">";
	  	  	os->width(14);
	  	  	*os << this->putLatticeVolume() / 3.0;
		}
  	}
  	else{
  		*os << "<VolumeOfUnitCell>";
  	  	os->width(14);
  	  	*os << this->putLatticeVolume();
  	}
  	*os << " </VolumeOfUnitCell>\n";

  	const SetOfFigureOfMerit& critical_value = this->putFiguresOfMerit();

	os->width(label_start); *os << "";
   	*os << "<FigureOfMeritWolff name=\"" << critical_value.putLabel_FigureOfMeritWolff() << "\">";
  	os->width(14);
   	*os << critical_value.putFigureOfMeritWolff();
	*os << " </FigureOfMeritWolff>\n";

	os->width(label_start);
   	*os << "" << "<FigureOfMeritWu name=\"" << critical_value.putLabel_FigureOfMeritWu() << "\">";
  	os->width(14);
   	*os << critical_value.putFigureOfMeritWu();
	*os << " </FigureOfMeritWu>\n";

	os->width(label_start);
   	*os << "" << "<ReversedFigureOfMeritWolff name=\"" << critical_value.putLabel_ReversedFigureOfMeritWolff() << "\">";
  	os->width(14);
   	*os << critical_value.putReversedFigureOfMerit();
	*os << " </ReversedFigureOfMeritWolff>\n";

	os->width(label_start);
   	*os << "" << "<SymmetricFigureOfMeritWolff name=\"" << critical_value.putLabel_SymmetricFigureOfMeritWolff() << "\">";
  	os->width(14);
   	*os << critical_value.putSymmetricFigureOfMerit();
	*os << " </SymmetricFigureOfMeritWolff>\n";

	os->width(label_start);
   	*os << "" << "<!-- Number of reflections up to the ";
   	*os << critical_value.putNumberOfReflectionsForFigureOfMerit() << "th reflection. (The multiplicity of peaks is considered.)-->\n";
	os->width(label_start);
   	*os << "" << "<NumberOfMillerIndicesInRange>";
  	os->width(14);
   	*os << critical_value.putContinuousNumberOfHKLInRange();
   	*os << " </NumberOfMillerIndicesInRange>\n";

   	os->width(label_start);
   	*os << "" << "<NumberOfIndexedPeaks>";
  	os->width(14);
   	*os << critical_value.putNumQobsAssociatedWithCloseHKL();
   	*os << " </NumberOfIndexedPeaks>\n";

   	os->width(label_start);
   	*os << "" << "<NumberOfPeaksNecessaryToResolveDominantZone>";
  	os->width(14);
   	*os << this->checkDominantZone();
   	*os << " </NumberOfPeaksNecessaryToResolveDominantZone>\n\n";



   	os->width(label_start);
	*os << "" << "<EquivalentLatticeCandidates>\n";
	*os << "" << "<!-- Almost equivalent unitcell parameters of different centring types. -->\n";
	label_start++;

	vector< pair< eBravaisType, pair< VecDat3<Double>, VecDat3<Double> > > > lattice_equiv;
	this->putEquivalentLatticeConstantsDegreeWithOtherCentring(abc_axis, rh_axis, resol, lattice_equiv);

		for(vector< pair< eBravaisType, pair< VecDat3<Double>, VecDat3<Double> > > >::const_iterator it=lattice_equiv.begin(); it<lattice_equiv.end(); it++)
		{
			os->width(label_start); *os << "";
	      	*os << "<LatticeCandidate>\n";
			label_start++;

			os->width(label_start); *os << "";
	      	*os << "<Centring>";
			os->width(17);
			if( it->first == Rhombohedral && rh_axis == Rho_Axis )
			{
				*os << "Rhombohedral(rhombohedral-axis)";
			}
			else *os << put_centring_name( BravaisType(it->first, abc_axis, rh_axis).enumCentringType() );
	      	*os << " </Centring>\n";

			os->width(label_start);
	      	*os << "" << "<LatticeParameters>";
	      	os->width(14);
	      	*os << it->second.first[0];
	      	os->width(14);
	   	   	*os << it->second.first[1];
	      	os->width(14);
	   	   	*os << it->second.first[2];
	      	os->width(14);
	   	   	*os << it->second.second[0];
	      	os->width(14);
	   	   	*os << it->second.second[1];
	      	os->width(14);
	   	   	*os << it->second.second[2];
	      	*os << " </LatticeParameters>\n";

	   	   	label_start--;
			os->width(label_start); *os << "";
	      	*os << "</LatticeCandidate>\n";
		}

	label_start--;
	os->width(label_start); *os << "";
	*os << "</EquivalentLatticeCandidates>\n\n";
}


void LatticeFigureOfMerit::putLatticeConstantsDegree(const BravaisType& brat, const SymMat<Double>& S0,
		const eABCaxis& axis1,
		const eRHaxis& axis2, VecDat3<Double>& length_axis, VecDat3<Double>& angle_axis)
{
	SymMat<Double> S = S0;
	if( brat.enumBravaisType() == Rhombohedral && axis2 != brat.enumRHaxis() )
	{
		if( axis2 == Hex_Axis ) // Rho -> Hex.
		{
			static const FracMat matrix_rho2hex = FInverse3( transpose(BravaisType::putTransformMatrixFromPrimitiveToRhomHex() ) );
			S = transform_sym_matrix(matrix_rho2hex.mat, S)/(matrix_rho2hex.denom*matrix_rho2hex.denom);
		}
		else // if( axis2 == RhoAxis ) // Hex -> Rho.
		{
			static const NRMat<Int4> matrix_hex2rho = transpose( BravaisType::putTransformMatrixFromPrimitiveToRhomHex() );
			S = transform_sym_matrix(matrix_hex2rho, S);
		}
	}
	else if( brat.enumBravaisType() == Monoclinic_B )
	{
		const NRMat<Int4> this2output = put_transform_matrix_monoclinic_b(brat.enumABCaxis(), axis1);
		S = transform_sym_matrix(this2output, S);
	}

	calLatticeConstant( S, length_axis, angle_axis );
}



Int4 LatticeFigureOfMerit::checkDominantZone() const
{
	const vector<QData> qdata = VCData::putPeakQData();
	if( qdata.empty() )
	{
		if( this->enumLaueGroup() == Ci ) return 6;
		else if( this->enumLaueGroup() == C2h_X ||  this->enumLaueGroup() == C2h_Y ||  this->enumLaueGroup() == C2h_Z ) return 4;
		else if( this->enumLaueGroup() == D2h ) return 3;
		else if( this->enumLaueGroup() == D4h_Z || this->enumLaueGroup() == D31d_rho || this->enumLaueGroup() == D6h ) return 2;
		else if( this->enumLaueGroup() == Oh ) return 1;
		assert(false);
	}

	const SymMat<Double> S_super = this->putSellingReducedForm();
	const Double max_q = max(
						  max( max( S_super(0,0), S_super(1,1) ), max( S_super(2,2), S_super(3,3) ) ),
						  max( max( S_super(0,0) + S_super(1,1) + S_super(0,1)*2.0,
								 	 S_super(0,0) + S_super(2,2) + S_super(0,2)*2.0 ),
								 	 S_super(1,1) + S_super(2,2) + S_super(1,2)*2.0 ) );

	return distance( qdata.begin(), lower_bound( qdata.begin(), qdata.end(), QData( qdata.begin()->q + max_q, 0.0 ) ) ) + 1;
}

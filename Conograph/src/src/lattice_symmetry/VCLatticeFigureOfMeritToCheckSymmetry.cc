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
#include "../utility_func/chToDouble.hh"
#include "../utility_lattice_reduction/put_Buerger_reduced_lattice.hh"
#include "OutputInfo.hh"
#include "VCLatticeFigureOfMeritToCheckSymmetry.hh"

const string VCLatticeFigureOfMeritToCheckSymmetry::CS_LABEL[NUM_LS] =
		{	"01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14" };

VCLatticeFigureOfMeritToCheckSymmetry::VCLatticeFigureOfMeritToCheckSymmetry(const Double& rhs)
	: m_label(-1), m_latfom(rhs),
	  m_S_red( SymMat43_VCData( SymMat<VCData>(3), NRMat<Int4>(4,3) ) )
{
}


VCLatticeFigureOfMeritToCheckSymmetry::VCLatticeFigureOfMeritToCheckSymmetry(const BravaisType& brat,
		const SymMat43_VCData& S,
		const ePeakShiftFunctionType& type,
		const Double& wave_length,
		const vector<ZParawError>& peak_shift_param_rad)
	: m_label(-1),
	  m_S_red( SymMat43_VCData( SymMat<VCData>(3), NRMat<Int4>(4,3) ) )
{
	this->setLatticeConstants43(brat, S);
	m_latfom.setPeakShiftParamRadian(type, wave_length, peak_shift_param_rad);
	m_num_lattice_found = 0;
}


#ifdef DEBUG
static bool checkInitialLatticeParameters(
		const BravaisType& brat,
		const SymMat43_VCData& S_red)
{
	const SymMat<Double> dbl_S_red( chToDouble(S_red.first) );

	if( brat.enumLaueGroup() == Ci && brat.enumCentringType() == Prim )
	{
		assert( dbl_S_red(2,2)*0.9999 < dbl_S_red(1,1) && dbl_S_red(1,1)*0.9999 < dbl_S_red(0,0)
				&& fabs( dbl_S_red(0,1) ) * 1.9999 < dbl_S_red(1,1)
				&& fabs( dbl_S_red(0,2) ) * 1.9999 < dbl_S_red(2,2)
				&& fabs( dbl_S_red(1,2) ) * 1.9999 < dbl_S_red(2,2) );
	}
	else if( brat.enumLaueGroup() == C2h_Y && brat.enumCentringType() == Prim )
	{
		assert( 0.0 <= dbl_S_red(0,2) && dbl_S_red(2,2)*0.9999 < dbl_S_red(0,0)
				&& fabs( dbl_S_red(0,2) ) * 1.9999 < dbl_S_red(2,2) && fabs( dbl_S_red(0,2) ) * 1.9999 < dbl_S_red(0,0) );
	}
	else if( brat.enumLaueGroup() == C2h_Z && brat.enumCentringType() == Prim )
	{
		assert( 0.0 <= dbl_S_red(0,1) && dbl_S_red(1,1)*0.9999 < dbl_S_red(0,0)
				&& fabs( dbl_S_red(0,1) ) * 1.9999 < dbl_S_red(0,0) && fabs( dbl_S_red(0,1) ) * 1.9999 < dbl_S_red(1,1) );
	}
	else if( brat.enumLaueGroup() == C2h_X && brat.enumCentringType() == Prim )
	{
		assert( 0.0 <= dbl_S_red(1,2) && dbl_S_red(2,2)*0.9999 < dbl_S_red(1,1)
				&& fabs( dbl_S_red(1,2) ) * 1.9999 < dbl_S_red(1,1) && fabs( dbl_S_red(1,2) ) * 1.9999 < dbl_S_red(2,2) );
	}
	else if( brat.enumLaueGroup() == C2h_Y && brat.enumCentringType() == BaseZ )
	{
		assert( 0.0 <= dbl_S_red(0,2) && fabs( dbl_S_red(0,2) ) * 1.9999 < dbl_S_red(2,2) && fabs( dbl_S_red(0,2) ) * 0.9999 < dbl_S_red(0,0) );
	}
	else if( brat.enumLaueGroup() == C2h_Z && brat.enumCentringType() == BaseX )
	{
		assert( 0.0 <= dbl_S_red(0,1) && fabs( dbl_S_red(0,1) ) * 1.9999 < dbl_S_red(0,0) && fabs( dbl_S_red(0,1) ) * 0.9999 < dbl_S_red(1,1) );
	}
	else if( brat.enumLaueGroup() == C2h_X && brat.enumCentringType() == BaseY )
	{
		assert( 0.0 <= dbl_S_red(1,2) && fabs( dbl_S_red(1,2) ) * 1.9999 < dbl_S_red(1,1) && fabs( dbl_S_red(1,2) ) * 0.9999  < dbl_S_red(2,2) );
	}
	else if( brat.enumLaueGroup() == D2h
			&& brat.enumCentringType() != BaseX
			&& brat.enumCentringType() != BaseY
			&& brat.enumCentringType() != BaseZ )
	{
		assert( dbl_S_red(2,2)*0.9999 < dbl_S_red(1,1) && dbl_S_red(1,1)*0.9999 < dbl_S_red(0,0) );
	}

	const SymMat<VCData> S_super = transform_sym_matrix(S_red.second, S_red.first);
	assert( S_super(0,1) <= VCData()
			&& S_super(0,2) <= VCData()
			&& S_super(0,3) <= VCData()
			&& S_super(1,2) <= VCData()
			&& S_super(1,3) <= VCData()
			&& S_super(2,3) <= VCData() );

	SymMat<VCData> S_red_cp = S_red.first;
	cal_average_crystal_system(brat.enumLaueGroup(), S_red_cp);
	assert( S_red.first(0,0).putVecCoef() == S_red_cp(0,0).putVecCoef() );
	assert( S_red.first(1,1).putVecCoef() == S_red_cp(1,1).putVecCoef() );
	assert( S_red.first(2,2).putVecCoef() == S_red_cp(2,2).putVecCoef() );
	assert( S_red.first(0,1).putVecCoef() == S_red_cp(0,1).putVecCoef() );
	assert( S_red.first(0,2).putVecCoef() == S_red_cp(0,2).putVecCoef() );
	assert( S_red.first(1,2).putVecCoef() == S_red_cp(1,2).putVecCoef() );

	return true;
}
#endif


void VCLatticeFigureOfMeritToCheckSymmetry::setLatticeConstants43(const BravaisType& brat, const SymMat43_VCData& S)
{
	m_S_red = S;
	assert( checkInitialLatticeParameters(brat, m_S_red) );

	m_latfom.setLatticeConstants43(brat, SymMat43_Double(chToDouble(m_S_red.first), m_S_red.second));
	m_num_lattice_found = 0;
}




bool VCLatticeFigureOfMeritToCheckSymmetry::checkIfLatticeIsMonoclinic(const ePointGroup& epg_new,
		const Double& cv2,
		map< SymMat<VCData>, NRMat<Int4> >& ans) const
{
	ans.clear();

	SymMat<VCData> ans0 = m_S_red.first;
	cal_average_crystal_system(C2h_X, ans0);

	SymMat<VCData> S_red(3);
	NRMat<Int4> trans_mat2;
	if( check_equiv_m(ans0, m_S_red.first, cv2 ) )
	{
		if( epg_new == C2h_X )
		{
			S_red = ans0;
			trans_mat2 = m_S_red.second;
			putBuergerReducedMonoclinicP<VCData>(1, 2, S_red, trans_mat2);
		}
		else if( epg_new == C2h_Y )
		{
			S_red = transform_sym_matrix(put_matrix_YXZ(), ans0);
			trans_mat2 = mprod(m_S_red.second, put_matrix_YXZ());
			putBuergerReducedMonoclinicP<VCData>(0, 2, S_red, trans_mat2);
		}
		else // if( epg_new == C2h_Z )
		{
			S_red = transform_sym_matrix(put_matrix_YZX(), ans0);
			trans_mat2 = mprod(m_S_red.second, put_matrix_ZXY());
			putBuergerReducedMonoclinicP<VCData>(0, 1, S_red, trans_mat2);
		}
		ans.insert( SymMat43_VCData( S_red, trans_mat2) );
	}

	ans0 = m_S_red.first;
	cal_average_crystal_system(C2h_Y, ans0);
	if( check_equiv_m(ans0, m_S_red.first, cv2 ) )
	{
		if( epg_new == C2h_X )
		{
			S_red = transform_sym_matrix(put_matrix_YXZ(), ans0);
			trans_mat2 = mprod(m_S_red.second, put_matrix_YXZ());
			putBuergerReducedMonoclinicP<VCData>(1, 2, S_red, trans_mat2);
		}
		else if( epg_new == C2h_Y )
		{
			S_red = ans0;
			trans_mat2 = m_S_red.second;
			putBuergerReducedMonoclinicP<VCData>(0, 2, S_red, trans_mat2);
		}
		else // if( epg_new == C2h_Z )
		{
			S_red = transform_sym_matrix(put_matrix_XZY(), ans0);
			trans_mat2 = mprod(m_S_red.second, put_matrix_XZY());
			putBuergerReducedMonoclinicP<VCData>(0, 1, S_red, trans_mat2);
		}
		ans.insert( SymMat43_VCData( S_red, trans_mat2) );
	}

	ans0 = m_S_red.first;
	cal_average_crystal_system(C2h_Z, ans0);
	if( check_equiv_m(ans0, m_S_red.first, cv2 ) )
	{
		if( epg_new == C2h_X )
		{
			S_red = transform_sym_matrix(put_matrix_ZXY(), ans0);
			trans_mat2 = mprod(m_S_red.second, put_matrix_YZX());
			putBuergerReducedMonoclinicP<VCData>(1, 2, S_red, trans_mat2);
		}
		else if( epg_new == C2h_Y )
		{
			S_red = transform_sym_matrix(put_matrix_XZY(), ans0);
			trans_mat2 = mprod(m_S_red.second, put_matrix_XZY());
			putBuergerReducedMonoclinicP<VCData>(0, 2, S_red, trans_mat2);
		}
		else // if( epg_new == C2h_Z )
		{
			S_red = ans0;
			trans_mat2 = m_S_red.second;
			putBuergerReducedMonoclinicP<VCData>(0, 1, S_red, trans_mat2);
		}
		ans.insert( SymMat43_VCData( S_red, trans_mat2) );
	}

	return !( ans.empty() );
}


bool VCLatticeFigureOfMeritToCheckSymmetry::checkIfLatticeIsOrthorhombic(const Double& cv2,
							map< SymMat<VCData>, NRMat<Int4> >& ans) const
{
	ans.clear();

	const BravaisType& brat = m_latfom.putBravaisType();

	SymMat<VCData> ans0 = m_S_red.first;
	cal_average_crystal_system(D2h, ans0);
	if( check_equiv_m(ans0, m_S_red.first, cv2 ) )
	{
		if( brat.enumCentringType() == BaseX )
		{
			if( ans0(1,1) < ans0(2,2) )
			{
				ans.insert( SymMat43_VCData( transform_sym_matrix(put_matrix_ZYX(), ans0), mprod( m_S_red.second, put_matrix_ZYX() ) ) );
			}
			else
			{
				ans.insert( SymMat43_VCData( transform_sym_matrix(put_matrix_YZX(), ans0), mprod( m_S_red.second, put_matrix_ZXY() ) ) );
			}
		}
		else if( brat.enumCentringType() == BaseY )
		{
			if( ans0(0,0) < ans0(2,2) )
			{
				ans.insert( SymMat43_VCData( transform_sym_matrix(put_matrix_ZXY(), ans0), mprod( m_S_red.second, put_matrix_YZX() ) ) );
			}
			else
			{
				ans.insert( SymMat43_VCData( transform_sym_matrix(put_matrix_XZY(), ans0), mprod( m_S_red.second, put_matrix_XZY() ) ) );
			}
		}
		else if( brat.enumCentringType() == BaseZ )
		{
			if( ans0(0,0) < ans0(1,1) )
			{
				ans.insert( SymMat43_VCData( transform_sym_matrix(put_matrix_YXZ(), ans0), mprod( m_S_red.second, put_matrix_YXZ() ) ) );
			}
			else
			{
				ans.insert( SymMat43_VCData( transform_sym_matrix(put_matrix_XYZ(), ans0), m_S_red.second ) );
			}
		}
		else
		{
			NRMat<Int4> trans_mat = m_S_red.second;
			putBuergerReducedOrthorhombic(brat.enumCentringType(), ans0, trans_mat);
			ans.insert( SymMat43_VCData(ans0, trans_mat ) );
		}
		return true;
	}
	return false;
}


bool VCLatticeFigureOfMeritToCheckSymmetry::checkIfLatticeIsTetragonal(const Double& cv2,
		map< SymMat<VCData>, NRMat<Int4> >& ans) const
{
	ans.clear();

	SymMat<VCData> ans0 = m_S_red.first;
	cal_average_crystal_system(D4h_X, ans0);
	if( check_equiv_m(ans0, m_S_red.first, cv2 ) )
	{
		ans.insert( SymMat43_VCData(
				transform_sym_matrix(put_matrix_YZX(), ans0), mprod( m_S_red.second, put_matrix_ZXY() ) ) );
	}

	ans0 = m_S_red.first;
	cal_average_crystal_system(D4h_Y, ans0);
	if( check_equiv_m(ans0, m_S_red.first, cv2 ) )
	{
		ans.insert( SymMat43_VCData(
				transform_sym_matrix(put_matrix_XZY(), ans0), mprod( m_S_red.second, put_matrix_XZY() ) ) );
	}

	ans0 = m_S_red.first;
	cal_average_crystal_system(D4h_Z, ans0);
	if( check_equiv_m(ans0, m_S_red.first, cv2 ) )
	{
		ans.insert( SymMat43_VCData(ans0, m_S_red.second ) );
	}

	return !( ans.empty() );
}




bool VCLatticeFigureOfMeritToCheckSymmetry::checkIfLatticeIsHexagonal(const ePointGroup& epg_new, const Double& cv2,
		map< SymMat<VCData>, NRMat<Int4> >& ans) const
{
	ans.clear();
	const BravaisType& brat = m_latfom.putBravaisType();

	SymMat43_VCData ans2(SymMat<VCData>(3), NRMat<Int4>(3,3));

	if( brat.enumLaueGroup() == C2h_X )
	{
		ans2.first = transform_sym_matrix(put_matrix_YZX(), m_S_red.first);
		ans2.second = mprod( m_S_red.second, put_matrix_ZXY() );
	}
	else if( brat.enumLaueGroup() == C2h_Y )
	{
		ans2.first = transform_sym_matrix(put_matrix_ZXY(), m_S_red.first);
		ans2.second = mprod( m_S_red.second, put_matrix_YZX() );
	}
	else // if( brat.enumLaueGroup() == C2h_Z )
	{
		ans2.first = transform_sym_matrix(put_matrix_XYZ(), m_S_red.first);
		ans2.second = m_S_red.second;
	}

	if( ans2.first(0,1) < VCData() )
	{
		ans2.first(0,1) *= -1;
		ans2.second[0][0] *= -1;
		ans2.second[1][0] *= -1;
		ans2.second[2][0] *= -1;
	}

	SymMat<VCData> ans0 = ans2.first;
	cal_average_crystal_system(epg_new, ans2.first);
	if( check_equiv_m(ans2.first, ans0, cv2 ) )
	{
		ans.insert( ans2 );
		return true;
	}
	else return false;
}


bool VCLatticeFigureOfMeritToCheckSymmetry::checkLatticeSymmetry(const ePointGroup& epg_new, const Double& cv2,
		map< SymMat<VCData>, NRMat<Int4> >& ans) const
{
	ans.clear();
	const BravaisType& brat = m_latfom.putBravaisType();
	if( epg_new == Ci || epg_new == brat.enumLaueGroup() )
	{
		ans.insert( m_S_red );
		return true;
	}

	if( epg_new == C2h_X || epg_new == C2h_Y ||  epg_new == C2h_Z )
	{
		assert( brat.enumLaueGroup() == Ci );
		assert( brat.enumCentringType() == Prim );

		return checkIfLatticeIsMonoclinic(epg_new, cv2, ans);
	}
	else if( epg_new == D4h_Z )
	{
		assert( brat.enumLaueGroup() == D2h );
		assert( brat.enumCentringType() == Prim
				|| brat.enumCentringType() == Inner );

		return checkIfLatticeIsTetragonal(cv2, ans);
	}
	else if( epg_new == D2h )
	{
		assert( brat.enumLaueGroup() != Ci || brat.enumCentringType() == Prim );
		assert( brat.enumLaueGroup() != C2h_Z || brat.enumCentringType() == BaseX );
		assert( brat.enumLaueGroup() != C2h_X || brat.enumCentringType() == BaseY );
		assert( brat.enumLaueGroup() != C2h_Y || brat.enumCentringType() == BaseZ );
		assert( brat.enumCentringType() != Rhom_hex );

		return checkIfLatticeIsOrthorhombic(cv2, ans);
	}
	else if( epg_new == D6h )
	{
		assert( brat.enumCentringType() == Prim );
		assert( brat.enumLaueGroup() == C2h_X
				|| brat.enumLaueGroup() == C2h_Y
				|| brat.enumLaueGroup() == C2h_Z );
		return checkIfLatticeIsHexagonal(epg_new, cv2, ans);
	}
	else
	{
		assert( epg_new == Oh );
		assert( brat.enumCentringType() == Prim
				|| brat.enumCentringType() == Inner
				|| brat.enumCentringType() == Face );

		SymMat43_VCData ans2 = m_S_red;
		cal_average_crystal_system(epg_new, ans2.first);
		if( check_equiv_m(ans2.first, m_S_red.first, cv2 ) )
		{
			ans.insert( ans2 );
			return true;
		}
	}
	return !(ans.empty());
}


void VCLatticeFigureOfMeritToCheckSymmetry::putLatticesOfHigherSymmetry(
		const ePointGroup& epg, const Double& cv2,
		vector<VCLatticeFigureOfMeritToCheckSymmetry>& lattice_result) const
{
	lattice_result.clear();
	map< SymMat<VCData>, NRMat<Int4> > S_red_tray;
	if( !this->checkLatticeSymmetry(epg, cv2, S_red_tray) ) return;

	const BravaisType& ebrat_original = this->putLatticeFigureOfMerit().putBravaisType();
	const eCentringType eblat = (ebrat_original.enumBravaisType()==Monoclinic_B?
									(epg==D31d_rho?Prim:(epg==D3d_1_hex?Rhom_hex:BaseZ)):ebrat_original.enumCentringType());

	const NRMat<Int4> matrix_min_to_sell = this->putInitialForm().second;

	SymMat<Double> S_super(4);
	NRMat<Int4> trans_mat(4,3);

	for(map< SymMat<VCData>, NRMat<Int4> >::const_iterator it=S_red_tray.begin(); it!=S_red_tray.end(); it++)
	{
		// S_super = it->second * it->first * Transpose(it->second) is close to
		// Delone-reduced form of the original lattice.
		S_super = transform_sym_matrix(it->second, chToDouble(it->first) );

		trans_mat = identity_matrix<Int4>(4);

		// S_super = trans_mat * it->second * it->first * Transpose(trans_mat * it->second).
		put_Selling_reduced_dim_less_than_4(S_super, trans_mat);
		moveSmallerDiagonalLeftUpper(S_super, trans_mat);

		lattice_result.push_back( VCLatticeFigureOfMeritToCheckSymmetry( BravaisType( pair<eCentringType, ePointGroup>(eblat, epg) ),
										SymMat43_VCData(it->first, mprod(trans_mat, it->second) ),
										this->putLatticeFigureOfMerit().putPeakShiftFunctionType(),
										this->putLatticeFigureOfMerit().putWaveLength(),
										this->putLatticeFigureOfMerit().putPeakShiftParamRadian() ) );
	}
}


void VCLatticeFigureOfMeritToCheckSymmetry::printLatticeInformation(
					const vector<VCLatticeFigureOfMeritToCheckSymmetry> lattice_result[],
					const OutputInfo outinfo[],
					const eABCaxis& abc_axis,
					const eRHaxis& rh_axis,
					const Double& resol,
					const Int4& label_start0,
					ostream* os) const
{
	m_latfom.printLatticeInformation(abc_axis, rh_axis, resol, label_start0, os);

  	Int4 label_start = label_start0;
	os->width(label_start);
	*os << "" << "<NumberOfLatticesInNeighborhood>";
  	os->width(14);
   	*os << this->putNumberOfLatticesInNeighborhood();
  	*os << " </NumberOfLatticesInNeighborhood>\n\n";
}

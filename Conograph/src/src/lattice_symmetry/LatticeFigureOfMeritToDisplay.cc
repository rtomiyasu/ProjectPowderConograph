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
#include "gather_q_of_3D_lattice.hh"
#include "../utility_func/zmath.hh"
#include "../zlog/zlog.hh"
#include "../ControlParam.hh"
#include "../PeakPosData.hh"
#include "LatticeFigureOfMeritToDisplay.hh"

LatticeFigureOfMeritToDisplay::LatticeFigureOfMeritToDisplay()
	: m_type_of_reflection_conditions(-1), m_showsTicks(false)
{
}


void LatticeFigureOfMeritToDisplay::setTypeOfSystematicAbsences(const Int4& arg)
{
	assert(arg < putNumberOfTypesOfSystematicAbsences(this->putLatticeFigureOfMerit().putBravaisType()));
	m_type_of_reflection_conditions = arg;
};


void LatticeFigureOfMeritToDisplay::putStandardMillerIndicesToFit(vector< VecDat3<Int4> >& hkl_to_fit,
vector<bool>& fix_or_fit_flag) const
{
	const Int4 num_Q = m_associated_hkl_tray.size();
	const Double cv2 = m_latfom.putCriticalValueSquare();

	fix_or_fit_flag.resize(num_Q);
	hkl_to_fit.resize(num_Q);

	vector< multimap<Double, vector<HKL_Q>::const_iterator> >::const_iterator it = m_associated_hkl_tray.begin();
	Double diff;
	for(Int4 i=0; i<num_Q; i++, it++)
	{
		if( it->empty() )
		{
			hkl_to_fit[i] = 0;
			fix_or_fit_flag[i] = false;
			continue;
		}

		diff = m_qdata[i].q - it->begin()->second->Q();
		hkl_to_fit[i] = it->begin()->second->HKL();

		if( diff * diff <= cv2 * m_qdata[i].q_var ) fix_or_fit_flag[i] = true;
		else fix_or_fit_flag[i] = false;
	}
}


void LatticeFigureOfMeritToDisplay::resetMillerIndicesInRange(const Int4& num_fit_data)
{
	assert( num_fit_data <= (Int4)m_qdata.size() );

	const Double cv2 = m_latfom.putCriticalValueSquare();

	range<Double> qobs_range;
	qobs_range.begin = 0.0;
	qobs_range.end = m_qdata.rbegin()->q + sqrt(cv2 * m_qdata.rbegin()->q_var);

	m_latfom.putMillerIndicesInRange(qobs_range.end, m_type_of_reflection_conditions, m_cal_hkl_tray);
	m_associated_hkl_tray.clear();
	m_associated_hkl_tray.resize(num_fit_data);

	vector< multimap<Double, vector<HKL_Q>::const_iterator> >::iterator it = m_associated_hkl_tray.begin();
	vector<HKL_Q>::const_iterator it2_begin, it2_end;
	for(Int4 i=0; i<num_fit_data; i++, it++)
	{
		const Double diff = sqrt(cv2 * m_qdata[i].q_var);
		it2_begin = lower_bound( m_cal_hkl_tray.begin(), m_cal_hkl_tray.end(), HKL_Q(0, m_qdata[i].q - diff) );
		it2_end = upper_bound( it2_begin, (vector<HKL_Q>::const_iterator)m_cal_hkl_tray.end(), HKL_Q(0, m_qdata[i].q + diff) );
		if( it2_begin < it2_end )
		{
			for(vector<HKL_Q>::const_iterator it2=it2_begin; it2<it2_end; it2++)
			{
				it->insert( multimap<Double, vector<HKL_Q>::const_iterator>::value_type(fabs(m_qdata[i].q-it2->Q()), it2) );
			}
		}
		else
		{
			if( m_cal_hkl_tray.begin() < it2_begin )
			{
				it->insert( multimap<Double, vector<HKL_Q>::const_iterator>::value_type(fabs(m_qdata[i].q - (it2_begin-1)->Q()), it2_begin-1) );
			}
			if( it2_begin < m_cal_hkl_tray.end() )
			{
				it->insert( multimap<Double, vector<HKL_Q>::const_iterator>::value_type(fabs(m_qdata[i].q - it2_begin->Q()), it2_begin) );
			}
		}
	}
}



Int4 LatticeFigureOfMeritToDisplay::fitLatticeParameter(const PeakPosData& pdata,
		const vector<etype_ID>& fitflag, const Int4& Max_ITNUM, const Double& limiter)
{
	vector< VecDat3<Int4> > hkl_to_fit = m_hkl_to_fit;
	vector<bool> fix_or_fit_flag = m_fix_or_fit_flag;	// 0:fix, 1:fit.
	Int4 itnum = 0;
	try{
		do{
			const Double old_fom = m_latfom.putFiguresOfMerit().putFigureOfMeritWolff();
ZLOG_INFO( "\n" + num2str( ++itnum ) + ") Initial " + m_latfom.putFiguresOfMerit().putLabel_FigureOfMeritWolff() + " = " + num2str(old_fom) + "\n" );

			vector<QData> qdata_modified;
			pair<bool, ZErrorMessage> ans = m_latfom.fitLatticeParameter(m_qdata,
																	hkl_to_fit, fix_or_fit_flag, true,
																	pdata, fitflag, Max_ITNUM,
																	VCData::putPeakQData().size(), qdata_modified);

			if( ans.second.putErrorType() != ZErrorNoError )
			{
				throw ans.second;
			}
			else if( !ans.first )
			{
ZLOG_INFO( num2str( m_latfom.putFiguresOfMerit().putLabel_FigureOfMeritWolff() ) + " was not improved.\n\n" );
				break;
			}
			else
			{
				// Optimization has succeeded.
				const Double& new_fom = m_latfom.putFiguresOfMerit().putFigureOfMeritWolff();
ZLOG_INFO( num2str( m_latfom.putFiguresOfMerit().putLabel_FigureOfMeritWolff() ) + " of optimized solution = " + num2str(new_fom) + "\n" );
				m_qdata = qdata_modified;
				resetMillerIndicesInRange(m_associated_hkl_tray.size());
				putStandardMillerIndicesToFit(hkl_to_fit, fix_or_fit_flag);
				if( (new_fom - old_fom) < old_fom * limiter ) break;
			}
		} while(true);
	}
	catch(const ZErrorMessage& zerr)
	{
		ZLOG_ERROR( zerr.printErrorLog() );
	}

	if( itnum > 2 )
	{
		m_hkl_to_fit = hkl_to_fit;
		m_fix_or_fit_flag = fix_or_fit_flag;
	}

	return itnum-1;
}


void LatticeFigureOfMeritToDisplay::putCalculatedPeakPosInRange(const ControlParam& cdata,
		Vec_DP& cal_pos_tray) const
{
	static const Double MAX_DP = numeric_limits<Double>::max();
	const vector<ZParawError>& peak_shift_param_rad = m_latfom.putPeakShiftParamRadian();

	cal_pos_tray.clear();
	if( cdata.IsAngleDispersion() )
	{
		for(vector<HKL_Q>::const_iterator it=m_cal_hkl_tray.begin(); it<m_cal_hkl_tray.end(); it++)
		{
			if( cdata.putWaveLength() * cdata.putWaveLength() * it->Q() > 4.0 ) break;
			cal_pos_tray.push_back( cal_theta2_deg(peak_shift_param_rad, cdata.putWaveLength(), sqrt(it->Q()) ) );
		}

		const Int4 isize = m_cal_hkl_tray.size();
		cal_pos_tray.resize(isize, -MAX_DP);
	}
	else
	{
		for(vector<HKL_Q>::const_iterator it=m_cal_hkl_tray.begin(); it<m_cal_hkl_tray.end(); it++)
		{
			cal_pos_tray.push_back( put_polynomial_value(cdata.putConvParam(), 1.0 / sqrt(it->Q()) ) );
		}
	}
}


void LatticeFigureOfMeritToDisplay::putPeakPosToFit(const ControlParam& cdata,
		Vec_DP& cal_q_tray, Vec_DP& cal_pos_tray) const
{
	const Vec_DP& conv_param = cdata.putConvParam();
	const vector<ZParawError>& peak_shift_param_rad = m_latfom.putPeakShiftParamRadian();
	const Double& wave_length = cdata.putWaveLength();

	const vector< VecDat3<Int4> >& hkl_to_fit = putMillerIndicesToFit();
	const size_t num_obsQ = hkl_to_fit.size();
	const SymMat<Double>& dbl_S = m_latfom.putOptimizedForm().first;

	cal_q_tray.clear();
	cal_pos_tray.clear();
	cal_q_tray.resize(num_obsQ);
	cal_pos_tray.resize(num_obsQ);
	for(size_t l=0; l<num_obsQ; l++)
   	{
   		const VecDat3<Int4>& hkl = hkl_to_fit.at(l);
   		const Double q_cal = norm(hkl, dbl_S);
   		cal_q_tray.at(l) = q_cal;

   		if( cdata.IsAngleDispersion() )
   		{
   			cal_pos_tray.at(l) = cal_theta2_deg(peak_shift_param_rad, wave_length, sqrt(q_cal));
   		}
   		else
   		{
   			cal_pos_tray.at(l) = put_polynomial_value(conv_param, 1.0 / sqrt(q_cal));
   		}
   	}
}




void LatticeFigureOfMeritToDisplay::printIndexingResult(const ControlParam& cdata,
					const PeakPosData& pdata,
					const Int4& label_start0,
					ostream* os) const
{
	static const Double DegRad =  180.0 / PI();

	Int4 label_start = label_start0;
	os->width(label_start);
	*os << "" << "<IndexingResults>\n";
	label_start++;

	const Vec_DP& posdata = pdata.putPeakPosXData();
	const Vec_DP& posdata_err = pdata.putPeakWidthData();

	const vector<Int4>& pos_qindex = VCData::putPeakQIndex();

	const vector< VecDat3<Int4> >& hkl_to_fit = putMillerIndicesToFit();
	const vector<bool>& fix_or_fit_flag = putFittingIDs();
	const size_t num_obsQ = hkl_to_fit.size();

	Vec_DP cal_q_tray, cal_pos_tray;
	this->putPeakPosToFit(cdata, cal_q_tray, cal_pos_tray);

	os->width(label_start);
	*os << "" << "<!-- q_obs, q_err, q_cal, peak_pos, peak_width, pos_cal, hkl, fix_or_fit.-->\n";
	for(size_t l=0; l<num_obsQ; l++)
   	{
   		os->width(15);
   		*os << m_qdata.at(l).q;

   		os->width(15);
   		*os << sqrt(m_qdata.at(l).q_var);

   		os->width(15);
   		*os << cal_q_tray.at(l);

   		os->width(15);
   		*os << posdata[pos_qindex.at(l)];

  		os->width(15);
   		*os << posdata_err[pos_qindex.at(l)];

   		os->width(15);
		*os << cal_pos_tray.at(l);

   		const VecDat3<Int4>& hkl = hkl_to_fit.at(l);

   		os->width(15);
   		*os << "[" + num2str<Int4>( hkl[0] ) + ","
   				+ num2str<Int4>( hkl[1] ) + ","
   				+ num2str<Int4>( hkl[2] ) + "]";

   		os->width(5);
   	   	if( fix_or_fit_flag.at(l) ) *os << "1";
   	   	else *os << "0";

   		*os << endl;
   	}

	label_start--;
	os->width(label_start);
	*os << "" << "</IndexingResults>\n\n";

	if( !cdata.IsAngleDispersion() ) return;

	label_start = label_start0;
	os->width(label_start);
	*os << "" << "<OptimizedZeroShitParameters>\n";
	label_start++;

	vector<ZParawError> peak_shift_param_rad = m_latfom.putPeakShiftParamRadian();
	const Int4 param_num = peak_shift_param_rad.size();
	for(Int4 i=0; i<param_num; i++)
	{
		os->width(label_start);
		*os << "" << "<Parameter>";
		os->width(15);
		*os << peak_shift_param_rad[i].value * DegRad;
		*os << "  </Parameter>";

		*os << "  <Error>";
		os->width(15);
		*os << peak_shift_param_rad[i].error * DegRad;
		*os << "  </Error>\n";
	}

	label_start--;
	os->width(label_start);
	*os << "" << "</OptimizedZeroShitParameters>\n";
}

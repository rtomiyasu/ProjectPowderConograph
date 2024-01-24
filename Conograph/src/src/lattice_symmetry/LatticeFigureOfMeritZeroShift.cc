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
#include "../model_function/PeakPosModel.hh"
#include "../PeakPosData.hh"
#include "../chToqValue.hh"
#include "../utility_lattice_reduction/put_Buerger_reduced_lattice.hh"
#include "LatticeFigureOfMeritZeroShift.hh"


LatticeFigureOfMeritZeroShift::LatticeFigureOfMeritZeroShift()
{
	m_etype_peak_shift = kPeakShiftFunction_Type0;
	m_WlengthX = 0.0;
}


LatticeFigureOfMeritZeroShift::LatticeFigureOfMeritZeroShift(const Double& rhs) : LatticeFigureOfMerit(rhs)
{
}


ZErrorMessage LatticeFigureOfMeritZeroShift::setPeakShiftParamRadian(const ePeakShiftFunctionType& type,
		const Double& wave_length,
		const vector<ZParawError>& peak_shift_param_rad)
{
	m_etype_peak_shift = type;
	m_WlengthX = wave_length;

	assert( m_etype_peak_shift != kPeakShiftFunction_Type0 || peak_shift_param_rad.size() == 0 );
	assert( m_etype_peak_shift != kPeakShiftFunction_Type1 || peak_shift_param_rad.size() == 1 );

	m_peak_shift_param_rad = peak_shift_param_rad;
	return ZErrorMessage();
}


ZErrorMessage LatticeFigureOfMeritZeroShift::setPeakShiftParamRadian(
		const vector<QData>& qdata,
		const ePeakShiftFunctionType& type,
		const Double& wave_length,
		const vector<ZParawError>& peak_shift_param_rad,
		const PeakPosData& pdata,
		const Int4& qdata_modified_size, vector<QData>& qdata_modified)
{
	ZErrorMessage zerr = this->setPeakShiftParamRadian(type, wave_length, peak_shift_param_rad);
	if( zerr.putErrorType() != ZErrorNoError ) return zerr;

	qdata_modified.resize(qdata_modified_size);
	copy( qdata.begin(), qdata.begin()+qdata_modified_size, qdata_modified.begin() );

	if( m_etype_peak_shift != kPeakShiftFunction_Type0 )
	{
		const Vec_DP& posdata = pdata.putPeakPosXData();
		const Vec_DP& posdata_err = pdata.putPeakWidthData();
		const vector<Int4>& peak_qindex = VCData::putPeakQIndex();
		assert( 0 <= qdata_modified_size && qdata_modified_size <= (Int4)peak_qindex.size() );

		for(Int4 i=0; i<qdata_modified_size; i++)
		{
			zerr = chXDataToQData_AD(posdata[peak_qindex[i]], posdata_err[peak_qindex[i]],
										peak_shift_param_rad[0].value, wave_length, qdata_modified[i]);
			if( zerr.putErrorType() != ZErrorNoError ) return zerr;
		}
	}

	return ZErrorMessage();
}


bool check_all_fixed(const vector<etype_ID>& peak_shift_fitflag)
{
	const Int4 num = peak_shift_fitflag.size();
	for(Int4 i=0; i<num; i++)
	{
		if( peak_shift_fitflag[i] != _ZRietveldIDFixed ) return false;
	}
	return true;
}


//  true : NormM has been improved.
//  false : NormM has not been improved.
pair<bool, ZErrorMessage> LatticeFigureOfMeritZeroShift::fitLatticeParameter(
		const vector<QData>& qdata,
		const vector< VecDat3<Int4> >& hkl_to_fit,
		const vector<bool>& fix_or_fit_flag,
		const bool& output_view_flag,
		const PeakPosData& pdata,
		const vector<etype_ID>& peak_shift_fitflag, const Int4& Max_ITNUM,
		const Int4& qdata_modified_size, vector<QData>& qdata_modified)
{
	assert( this->putFiguresOfMerit().putNumberOfReflectionsForFigureOfMerit() <= qdata_modified_size );
	assert( qdata_modified_size <= (Int4)qdata.size() );
	qdata_modified.resize(qdata_modified_size);
	copy( qdata.begin(), qdata.begin()+qdata_modified_size, qdata_modified.begin() );

	if( m_etype_peak_shift == kPeakShiftFunction_Type0
			|| check_all_fixed(peak_shift_fitflag) )
	{
		return LatticeFigureOfMerit::fitLatticeParameterLinear(qdata, hkl_to_fit, fix_or_fit_flag, output_view_flag);
	}

	static const Double RadDeg = PI() / 180.0; // = pi / 180.0.
	static const Double FWHM2Err = RadDeg / sqrt( 8.0 * log(2.0) );
	const Vec_DP& posdata = pdata.putPeakPosXData();
	const Vec_DP& posdata_fwhm = pdata.putPeakWidthData();
	const vector<Int4>& pos_qindex = VCData::putPeakQIndex();

	assert( hkl_to_fit.size() == fix_or_fit_flag.size() );
	assert( hkl_to_fit.size() <= posdata.size() );

	const Int4 isize = hkl_to_fit.size();
	
	Vec_DP ydata(isize), ydata_err(isize);
	Vec_BOOL nxfit(isize);
	Int4 data_num=0;
	
	for(Int4 i=0; i<isize; i++)
	{
		const Int4& peak_index = pos_qindex[i];
		ydata[i] = posdata[peak_index] * RadDeg;
		ydata_err[i] = posdata_fwhm[peak_index] * FWHM2Err;
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
	
	// Set initial parameters.
	const Int4 peak_shift_param_num = m_peak_shift_param_rad.size();
	const Int4 param_num = 6 + peak_shift_param_num;
	vector<Double> init_param(param_num);
	const SymMat<Double>& S_val = this->putOptimizedForm().first;
	init_param[0] = S_val(0,0);
	init_param[1] = S_val(1,1);
	init_param[2] = S_val(2,2);
	init_param[3] = S_val(1,2);
	init_param[4] = S_val(0,2);
	init_param[5] = S_val(0,1);

	LaueGroup lg(this->enumLaueGroup());
	Mat_DP_constr cmat;
	lg->putLatticeConstantFlag(cmat);

	cmat.resize(param_num);
	for(Int4 i=0,i2=6; i<peak_shift_param_num; i++,i2++)
	{
		init_param[i2] = m_peak_shift_param_rad[i].value;
		cmat[i2].ID = peak_shift_fitflag[i];
	}

	if( data_num <= countNumberOfIndependentParam(cmat.begin(),cmat.end()) )
	{
		return LatticeFigureOfMerit::fitLatticeParameterLinear(qdata, hkl_to_fit, fix_or_fit_flag, output_view_flag);
	}
	setIndex(cmat);

	PeakPosModel posModel(m_etype_peak_shift, m_WlengthX);
	posModel.setConstraint(&cmat[0]);
	Double chisq_all;
	pair<bool, ZErrorMessage> ans = posModel.setFittedParam(hkl_to_fit, ydata, ydata_err, nxfit,
										output_view_flag, 1.0e-3, Max_ITNUM, init_param, chisq_all);
	if( !(ans.first) )
	{
		return ans;
	}

	LatticeFigureOfMeritZeroShift new_lat(*this);
	SymMat<Double> S_red_optimized(3);
	vector<ZParawError> peak_shift_param_rad;
	posModel.putResult(S_red_optimized, peak_shift_param_rad);
	new_lat.setLatticeConstants(this->putBravaisType(), S_red_optimized);
	ans.second = new_lat.setPeakShiftParamRadian(qdata, m_etype_peak_shift, m_WlengthX, peak_shift_param_rad, pdata,
													qdata_modified_size, qdata_modified);
	if( ans.second.putErrorType() != ZErrorNoError )
	{
		return ans;
	}

	new_lat.setFigureOfMerit(this->putFiguresOfMerit().putNumberOfReflectionsForFigureOfMerit(), qdata_modified);

	if( cmpFOMdeWolff(new_lat, *this) )
	{
		*this = new_lat;
		return pair<bool, ZErrorMessage>(true, ZErrorMessage());
	}
	return pair<bool, ZErrorMessage>(false, ZErrorMessage());
}

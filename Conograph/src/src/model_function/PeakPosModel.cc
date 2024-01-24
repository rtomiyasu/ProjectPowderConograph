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
#include "../ControlParam.hh"
#include "../levenberg_marquardt/LemarqMethod.hh"
#include "profile_function/global_function/PeakShiftFunc.hh"
#include "PeakPosModel.hh"

static IPeakShiftFunc* generatePeakShiftFunc(const ePeakShiftFunctionType& type)
{
	assert( type != kPeakShiftFunction_Type0 );	// In this case, it is supposed to use linear optimization.

	if( type == kPeakShiftFunction_Type1 ) return new PeakShiftFunc_Type1();	// peak-shift is given by Bragg formula.
	else return new PeakShiftFunc_Type0();	// peak-shift equals zero.
}

PeakPosModel::PeakPosModel(const ePeakShiftFunctionType& type, const Double& wave_length)
	:  m_peak_shift_model( generatePeakShiftFunc(type) ),
		ParaNum(m_lat_model.putParamNum() + m_peak_shift_model->putParamNum())
{
	m_num_indep = m_lat_model.putNumberOfIndependentParam()
					+ m_peak_shift_model->putNumberOfIndependentParam();
	m_conv_model.setWaveLengthAD(wave_length);
}

PeakPosModel::~PeakPosModel()
{
	delete m_peak_shift_model;
}

Int4 PeakPosModel::putParamNum() const
{ 
	return ParaNum;
};

void PeakPosModel::putResult(SymMat<Double>& ans, vector<ZParawError>& param_peak_shift_rad) const
{
	m_lat_model.putResult(ans);

	m_peak_shift_model->putResult(param_peak_shift_rad);
}


Int4 PeakPosModel::putNumberOfIndependentParam() const
{
	return m_num_indep;
}


void PeakPosModel::putConstraint(constr_DP* const tray) const
{
	m_lat_model.putConstraint(tray);
	m_peak_shift_model->putConstraint(tray+m_lat_model.putParamNum());
}



void PeakPosModel::setConstraint(const constr_DP* tray)
{
	m_lat_model.setConstraint(tray);
	m_peak_shift_model->setConstraint(tray+m_lat_model.putParamNum());
	m_num_indep = m_lat_model.putNumberOfIndependentParam()
					+ m_peak_shift_model->putNumberOfIndependentParam();
}


void PeakPosModel::setCovariantMatrixAll(const SymMat<Double>& cov)
{
	m_lat_model.setCovariantMatrixAll(cov);
	SymMat<Double> cov2(m_peak_shift_model->putParamNum());
	copyPartOfCovariantMatrix(cov, m_lat_model.putParamNum(), putParamNum(), cov2);
	m_peak_shift_model->setCovariantMatrixAll(cov2);
}


Int4 PeakPosModel::setParam(Double* param, ostringstream& outPutString)
{
	Int4 ans = m_lat_model.setParamErrOut(param, outPutString);
	if( ans <= 0 ) return 0;
	Int4 ans2 = m_peak_shift_model->setParamErrOut(param+m_lat_model.putParamNum(), outPutString);
	if( ans2 <= 0 ) return 0;
	if( ans2 >= 2 ) return 2;
	return ans;
}


bool PeakPosModel::putFittingCoef(const VecDat3<Int4>& hkl, Double& lclparam, Double* dlclparam_dother, ostringstream& outPutString) const
{
	Double invdsq;
	if( !m_lat_model.putFittingCoefErrOut(hkl, invdsq, &dlclparam_dother[0], outPutString) ) return false;

	const Double dwidth = 1.0/sqrt(invdsq);
	Double theta2_rad;
	Double dtheta2_dinvdsq;
	Vec_DP dtheta2_dother;	// empty array.
	if( !m_conv_model.putFittingCoefErrOut(dwidth, theta2_rad, dtheta2_dinvdsq, dtheta2_dother, outPutString) ) return false;

	Double delta;
	Double ddelta_dtheta2_rad;
	if( !m_peak_shift_model->putFittingCoefErrOut(theta2_rad, delta, ddelta_dtheta2_rad,
						&dlclparam_dother[m_lat_model.putNumberOfIndependentParam()], outPutString) ) return false;

	lclparam = theta2_rad + delta;
	const Double dlclparam_dinvdsq = (1.0 + ddelta_dtheta2_rad)*dtheta2_dinvdsq;
	for(Int4 i=0; i<m_lat_model.putNumberOfIndependentParam(); i++)
	{
		dlclparam_dother[i] *= dlclparam_dinvdsq;
	}

	return true;
}


// If wt != 2, lclsterr is used as weights in least square method. 
// If wt = 2, lclsterr is not used,  all the weights equal 1. 
pair<bool, ZErrorMessage> PeakPosModel::setFittedParam(const vector< VecDat3<Int4> >& hkl,
		const vector<Double>& lclparam, const vector<Double>& lclvar,
		const vector<bool>& nxfit,
		const bool& output_view_flag, const Double& judge_conv,  const Int4& max_itnum,
		const vector<Double>& glbinit, Double& chisq_all)
{
	Int4 itnum;
	
	LemarqMethod marq;
	marq.setParam(output_view_flag, judge_conv, max_itnum, 0.001);
	
	MarquardtFmodelBase< VecDat3<Int4> > LMFunc(*this);
	ZErrorMessage zerr = LMFunc.addHistogram(hkl, lclparam, lclvar, nxfit);
	if( zerr.putErrorType() != ZErrorNoError ) return pair<bool, ZErrorMessage>(false, zerr);
	
	Vec_DP chisq;
	Int4 reason_terminate_LM;
	return marq.execute< VecDat3<Int4> >(glbinit, LMFunc, itnum, chisq_all, chisq, reason_terminate_LM);
}

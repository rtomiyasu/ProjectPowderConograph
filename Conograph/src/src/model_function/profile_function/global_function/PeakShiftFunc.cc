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
#include<cmath>
#include "../../../utility_func/zmath.hh"
#include "PeakShiftFunc.hh"

// delta 2theta = 0;
PeakShiftFunc_Type0::PeakShiftFunc_Type0()
{
}

PeakShiftFunc_Type0::~PeakShiftFunc_Type0()
{
}


void PeakShiftFunc_Type0::putResult(vector<ZParawError>& param) const
{
	param.clear();
}


Int4 PeakShiftFunc_Type0::putParamNum() const
{
	return ParaNum;
}

void PeakShiftFunc_Type0::setCovariantMatrixAll(const SymMat<Double>& cov)
{
}


// set peakshift parameters.
Int4 PeakShiftFunc_Type0::setParam(Double* param, ostringstream& outPutString)
{
	return 1;
}


// theta2 : 2theta (deg)
// dlclparam_dother : d(2theta)/d(Z), d(2theta)/d(Ds), d(2theta)/d(Ts)
bool PeakShiftFunc_Type0::putFittingCoef(const Double& theta2, Double& delta, Double* ddelta_dother, ostringstream& outPutString) const
{
	delta = 0.0;

	return true;
}


// Calculate d(2theta)/dS00, d(2theta)/dS11, ....
bool PeakShiftFunc_Type0::putFittingCoef(const Double& theta2,
Double& delta, Double& ddelta_dtheta2, Double* ddelta_dother, ostringstream& outPutString) const
{
	// peak shift ( radian )
	delta = 0.0;
	ddelta_dtheta2 = 0.0;
	return true;
}




// delta 2theta = Z.
PeakShiftFunc_Type1::PeakShiftFunc_Type1()
{
	m_param.resize(ParaNum, 0.0);
	m_fitflag.resize(ParaNum, _ZRietveldIDFixed);
}

PeakShiftFunc_Type1::~PeakShiftFunc_Type1()
{
}


void PeakShiftFunc_Type1::putResult(vector<ZParawError>& param) const
{
	param = m_param;
}


Int4 PeakShiftFunc_Type1::putParamNum() const
{
	return ParaNum;
}

void PeakShiftFunc_Type1::setCovariantMatrixAll(const SymMat<Double>& cov)
{
	for(Int4 k=0; k<ParaNum; k++) m_param[k].error = sqrt_d( cov(k,k) );
}


// set peakshift parameters Z.
Int4 PeakShiftFunc_Type1::setParam(Double* param, ostringstream& outPutString)
{
	// Z.
	m_param[0] = param[0];
	return 1;
}


// theta2 : 2theta (deg)
// dlclparam_dother : d(2theta)/d(Z).
bool PeakShiftFunc_Type1::putFittingCoef(const Double& theta2_rad, Double& delta, Double* ddelta_dother, ostringstream& outPutString) const
{
	// peak shift ( radian )
	delta = -m_param[0].value;

	Int4 j=0;
	if( m_fitflag[0] != _ZRietveldIDFixed ) ddelta_dother[j++] = -1.0;

	return true;
}


bool PeakShiftFunc_Type1::putFittingCoef(const Double& theta2_rad,
Double& delta, Double& ddelta_dtheta2_rad, Double* ddelta_dother, ostringstream& outPutString) const
{
	// peak shift ( radian )
	delta = -m_param[0].value;
	ddelta_dtheta2_rad = 0.0;

	Int4 j=0;
	if( m_fitflag[0] != _ZRietveldIDFixed ) ddelta_dother[j++] = -1.0;

	return true;
}

void PeakShiftFunc_Type1::setConstraint(const constr_DP* fitflag)
{
	for(Int4 j=0; j<ParaNum; j++)
	{
		if( fitflag[j].ID == _ZRietveldIDFixed )
		{
			m_fitflag[j] = _ZRietveldIDFixed;
		}
		else  m_fitflag[j] = _ZRietveldIDVary;
	}
}

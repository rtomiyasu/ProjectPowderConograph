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
#include "GlbBraggDiffract.hh"

// f(d) is an equation of diffraction.
// Calculate 2theta with d and lambda.
// 2theta = 2*arcsin(lambda/2*d)
GlbBraggDiffract::GlbBraggDiffract()
{
	m_wlambda = 1.0;
}


GlbBraggDiffract::~GlbBraggDiffract()
{
}


void GlbBraggDiffract::putResult(vector<ZParawError>& param) const
{
	param.clear();
}

Int4 GlbBraggDiffract::putParamNum() const
{
	return ParaNum;
}

void GlbBraggDiffract::setCovariantMatrixAll(const SymMat<Double>& cov)
{
//	for(Int4 k=0; k<ParaNum; k++) m_param[k].error = sqrt_d( cov(k,k) );
}


void GlbBraggDiffract::setWaveLengthAD(const Double& wlambda)
{
	m_wlambda = wlambda;
}


Int4 GlbBraggDiffract::setParam(Double* param, ostringstream& outPutString)
{
	return 1;
}

// theta2 : 2theta (rad)
bool GlbBraggDiffract::putFittingCoef(const Double& dwidth, Double& theta2_rad, Double* dtheta2_dother, ostringstream& outPutString) const
{
	Double t = m_wlambda / (2.0 * dwidth);
	if( t > 1.0 || t < -1.0 ) return false;
	theta2_rad = 2.0 * asin( t );
	return true;
}


// Calculate d(2theta)/dS00, d(2theta)/dS11, ....
bool GlbBraggDiffract::putFittingCoef(const Double& dwidth,
		Double& theta2_rad, Double& dtheta2_dinvdsq, Vec_DP& dtheta2_dother, ostringstream& outPutString) const
{
	dtheta2_dother.clear();

	const Double t = m_wlambda / (2.0 * dwidth);
	if( t > 1.0 || t < -1.0 ) return false;
	theta2_rad = 2.0 * asin( t );

	// d(2theta) / dwidth = 2.0 / sqrt(1 - t^{2}) * t * (-1 / dwidth)
	//                    = -2.0 / (sqrt(1/(t*t) - 1) * dwidth).
	// d(2theta) / d(1/d^2) = -2.0 / (sqrt(1/(t*t) - 1) * dwidth) / (-2.0/dwidth^3).
	dtheta2_dinvdsq = dwidth * dwidth / sqrt(1.0/(t*t) - 1.0);

	return true;
}

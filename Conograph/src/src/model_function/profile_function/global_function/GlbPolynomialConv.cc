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
#include "GlbPolynomialConv.hh"
#include "../../../utility_func/covar_matrix.hh"

// f(d) is a polynomial of order Polyorder.
// m_param[1] >= 0
GlbPolynomialConv::GlbPolynomialConv(const Int4& order) : PolyOrder(order), ParaNum(order+1)
{
	assert(PolyOrder > 0);

	m_param.resize(ParaNum, 0.0);
	m_dparam.resize(PolyOrder, 0.0);
}

GlbPolynomialConv::~GlbPolynomialConv()
{
}

void GlbPolynomialConv::putResult(vector<ZParawError>& param) const
{
	param = m_param;
}

Int4 GlbPolynomialConv::putParamNum() const
{
	return ParaNum;
}

void GlbPolynomialConv::setCovariantMatrixAll(const SymMat<Double>& cov)
{
	for(Int4 k=0; k<ParaNum; k++) m_param[k].error = sqrt_d( cov(k,k) );
}

inline bool GlbPolynomialConv::checkParam(Double* param)
{
	if(param[1] < 0.0) param[1] = 0.0;
	return true;
}

Int4 GlbPolynomialConv::setParam(Double* param)
{
	if( !checkParam(param) ) return 0;
	
	for(Int4 k=0; k<ParaNum; k++) m_param[k] = param[k];
	for(Int4 k=1; k<ParaNum; k++) m_dparam[k-1] = k*param[k];

	return 1;
}

bool GlbPolynomialConv::putFittingCoef(const Double& dwidth, Double& lclparam, Double* dlclparam_dother) const
{
	Int4 k = ParaNum-1;
	lclparam = m_param[k].value;
	for(; k>0;)
	{
		--k;
		lclparam = lclparam*dwidth + m_param[k].value;
	}

	dlclparam_dother[0]=1.0;
	for(Int4 k=1; k<=PolyOrder; k++) dlclparam_dother[k]=dlclparam_dother[k-1]*dwidth;
	return true;
}

bool GlbPolynomialConv::putFittingCoef(const Double& dwidth,
Double& lclparam, Double& dlclparam_dinvdsq, Vec_DP& dlclparam_dother) const
{
	dlclparam_dother.resize(PolyOrder+1);
	if( !putFittingCoef(dwidth, lclparam, &dlclparam_dother[0]) ) return false;

	Double dlclparam_ddwidth = put_polynomial_value(m_dparam, dwidth);
	Double ddwidth_dinvdsq = -0.5*dwidth*dwidth*dwidth;

	dlclparam_dinvdsq = dlclparam_ddwidth * ddwidth_dinvdsq;
	return true;
}

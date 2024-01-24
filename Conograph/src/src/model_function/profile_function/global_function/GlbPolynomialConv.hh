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
#ifndef GLBPOLYNOMIALCONV_HH_
#define GLBPOLYNOMIALCONV_HH_

#include "../../../RietveldAnalysisTypes.hh"
#include "IConvFunc.hh"

using namespace std;

class GlbPolynomialConv : public IConvFunc
{
private:
	const Int4 PolyOrder; 
	const Int4 ParaNum; 
	vector<ZParawError> m_param;
	Vec_DP m_dparam;

	inline bool checkParam(Double* param);

protected:
	// Sets the values of parameters by the argument.
    Int4 setParam(Double*);
	bool putFittingCoef(const Double&, Double&, Double*) const;
	bool putFittingCoef(const Double&, Double&, Double&, Vec_DP&) const;

public:
	GlbPolynomialConv(const Int4&);
	virtual ~GlbPolynomialConv();

	Int4 putParamNum() const;
	void putResult(vector<ZParawError>&) const;
	void setCovariantMatrixAll(const SymMat<Double>&); 
};

#endif /*GLBPOLYNOMIALCONV_HH_*/

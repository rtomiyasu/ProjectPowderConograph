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
#ifndef GLBBRAGGDIFFRACT_HH_
#define GLBBRAGGDIFFRACT_HH_

#include "../../../RietveldAnalysisTypes.hh"
#include "IConvFunc.hh"

// f(d) = arcsin( lambda / 2*dwidth );
class GlbBraggDiffract : public IConvFunc
{
private:
	static const Int4 ParaNum = 0;
	Double m_wlambda; // wavelength

protected:
	// Sets the values of parameters by the argument.
    Int4 setParam(Double*, ostringstream& outPutString);

	// Return the value(Ymod) of the model function at X.
	// If dYdX is not NULL, also return the derivatives at X by dYdX.
	bool putFittingCoef(const Double& X, Double& Ymod, Double* dYdX, ostringstream& outPutString) const;
	bool putFittingCoef(const Double& dwidth,
							Double& lclparam, Double& dlclparam_dinvdsq, Vec_DP& dlclparam_dother, ostringstream& outPutString) const;

public:
	GlbBraggDiffract();
	virtual ~GlbBraggDiffract();

	Int4 putParamNum() const;
	void putResult(vector<ZParawError>&) const;
	void setWaveLengthAD(const Double&);
	void setCovariantMatrixAll(const SymMat<Double>&); 
};

#endif /*GLBBRAGGDIFFRACT_HH_*/

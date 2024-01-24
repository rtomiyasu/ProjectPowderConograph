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
#ifndef _IConvFunc_hh_
#define _IConvFunc_hh_

#include "IGlobalFunc.hh"

// Interface class of global model functions for PhaseModelFunction class.
class IConvFunc : public IGlobalFunc
{
private:
	virtual bool putFittingCoef(const Double& X,
									Double& Ymod, 
									Double& dY_dinvdsq, Vec_DP& dYdX, ostringstream& outPutString) const = 0;

public:
    virtual ~IConvFunc(){};

	// Return the number of all the parameters.(This includes the parameters fixed.)
	virtual Int4 putParamNum() const = 0;

	// Return the value(Ymod) of the model function at X.
	// Also return the derivatives by S using the Miller index hkl.
	// ( S is the quadratic form determined by the co-lattice constants. )
	// Unless this function is not newly overloaded, dYdS equal 0 for the calculation speed.

    // In this class, the flag on fix/fit is supposed to be always true; 
	// Returns the flags on fix/fit.
	void putFitFlag(vector<etype_ID>& fitflag) const { fitflag.clear(); fitflag.resize(putParamNum(), _ZRietveldIDVary); };
	// Set the constraints which indicates which parameters are fixed or independently fit.
	void setConstraint(const constr_DP*){ return; };

	inline bool putFittingCoefErrOut(const Double& X, Double& Ymod, Double& dY_dinvdsq, Vec_DP& dYdX, ostringstream& outPutString) const;
};


inline bool IConvFunc::putFittingCoefErrOut(const Double& X, Double& Ymod, Double& dY_dinvdsq, Vec_DP& dYdX, ostringstream& outPutString) const
{
	if( putFittingCoef(X, Ymod, dY_dinvdsq, dYdX, outPutString) ) return true;
	outPutString << string(typeid(*this).name()) + "::putFittingCoef failed.\n";

	return false;
}

#endif

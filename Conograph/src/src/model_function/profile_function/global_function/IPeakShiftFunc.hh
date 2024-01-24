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
#ifndef IPEAKSHIFTFUNC_HH_
#define IPEAKSHIFTFUNC_HH_

#include <typeinfo>
#include "../../../RietveldAnalysisTypes.hh"
#include "../../../zparam/ZParawError.hh"
#include "IGlobalFunc.hh"


// Class of Peak shift(delta 2theta);
class IPeakShiftFunc : public IGlobalFunc
{
protected:
	virtual bool putFittingCoef(const Double& theta2, Double& delta, Double& ddelta_dtheta2,
										Double* ddelta_dparam, ostringstream& outPutString) const = 0;

public:
	virtual ~IPeakShiftFunc(){};

	inline bool putFittingCoefErrOut(const Double& theta2, Double& delta, Double& ddelta_dtheta2,
										Double* ddelta_dparam, ostringstream& outPutString) const;
};


inline bool IPeakShiftFunc::putFittingCoefErrOut(const Double& X, Double& Ymod, Double& dY_dinvdsq, Double* dYdX, ostringstream& outPutString) const
{
	if( putFittingCoef(X, Ymod, dY_dinvdsq, dYdX, outPutString) ) return true;
	outPutString << string(typeid(*this).name()) + "::putFittingCoef failed.\n";

	return false;
}

#endif /*IPEAKSHIFTFUNC_HH_*/

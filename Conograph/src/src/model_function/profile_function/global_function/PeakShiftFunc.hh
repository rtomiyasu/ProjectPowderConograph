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
#ifndef PEAKSHIFTFUNC_HH_
#define PEAKSHIFTFUNC_HH_

#include "../../../RietveldAnalysisTypes.hh"
#include "IPeakShiftFunc.hh"

// delta 2theta = 0.
class PeakShiftFunc_Type0 : public IPeakShiftFunc
{
private:
	static const Int4 ParaNum = 0;
	
public:
	PeakShiftFunc_Type0();
	virtual ~PeakShiftFunc_Type0();

	Int4 putParamNum() const;
	void putResult(vector<ZParawError>& peak_shift_rad) const;
	Int4 setParam(Double* param, ostringstream& outPutString);
	void setCovariantMatrixAll(const SymMat<Double>&); 
	bool putFittingCoef(const Double&, Double&, Double*, ostringstream& outPutString) const;
	bool putFittingCoef(const Double&, Double&, Double&, Double* dlclparam_dother, ostringstream& outPutString) const;

	// Set the constraints which indicates which parameters are fixed or independently fit.
	void setConstraint(const constr_DP*){};
	void putFitFlag(vector<etype_ID>& fitflag) const { fitflag.clear(); };
};



// delta 2theta = Z;
class PeakShiftFunc_Type1 : public IPeakShiftFunc
{
private:
	static const Int4 ParaNum = 1;

	vector<ZParawError> m_param;
	vector<etype_ID> m_fitflag;
	
public:
	PeakShiftFunc_Type1();
	virtual ~PeakShiftFunc_Type1();

	Int4 putParamNum() const;
	void putResult(vector<ZParawError>&) const;
	Int4 setParam(Double* param, ostringstream& outPutString);
	void setCovariantMatrixAll(const SymMat<Double>&); 
	bool putFittingCoef(const Double&, Double&, Double*, ostringstream& outPutString) const;
	bool putFittingCoef(const Double&, Double&, Double&, Double* dlclparam_dother, ostringstream& outPutString) const;

	// Set the constraints which indicates which parameters are fixed or independently fit.
	void setConstraint(const constr_DP*);
	void putFitFlag(vector<etype_ID>& fitflag) const { fitflag = m_fitflag; };
};

#endif /*PeakShiftFunc_HH_*/

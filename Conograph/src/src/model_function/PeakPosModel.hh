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
#ifndef _PeakPosModel_hh_
#define _PeakPosModel_hh_

// PeakPosModelModel.hh

#include "LatticeDistanceModel.hh"
#include "profile_function/global_function/enumPeakShiftFunctionType.hh"
#include "profile_function/global_function/GlbBraggDiffract.hh"
#include "profile_function/global_function/IPeakShiftFunc.hh"

class ControlParam;
class GlbBraggDiffract;

// Interface class for LemarqMethod class.
class PeakPosModel : public IMarquardtFmodel< VecDat3<Int4> >
{
private:
	LatticeDistanceModel m_lat_model;
	GlbBraggDiffract m_conv_model;
	IPeakShiftFunc* m_peak_shift_model;

	const Int4 ParaNum;

	Int4 m_num_indep;

protected:
    Int4 setParam(Double*, ostringstream& outPutString);

	// Return the value(Ymod) of the model function at X.
	// If dYdX is not NULL, the entries of dYdP are the derivatives at X by the fitting parameters.
	bool putFittingCoef(const VecDat3<Int4>& X, Double& Ymod, Double* dYdP, ostringstream& outPutString) const;

	
public:
	PeakPosModel(const ePeakShiftFunctionType& type, const Double& wave_length);
	virtual ~PeakPosModel();
	
	Int4 putParamNum() const;
	
	void putResult(SymMat<Double>& ans, vector<ZParawError>& peak_shifts) const;
	
	void setCovariantMatrixAll(const SymMat<Double>&); 

    // Return the number of parameters independently fit.
    Int4 putNumberOfIndependentParam() const;

	// Returns constraints on the parameters.
	void putConstraint(constr_DP* const) const;

	// Set the constraints which indicates which parameters are fixed or independently fit.
	void setConstraint(const constr_DP*);

	pair<bool, ZErrorMessage> setFittedParam(const vector< VecDat3<Int4> >& hkl,
						const vector<Double>& lclparam,
						const vector<Double>& lclvar,
						const vector<bool>& nxfit,
						const bool& output_view_flag,
						const Double&,  const Int4&, const vector<Double>&, Double& chisq_all);

	LatticeDistanceModel& putLatticeDistanceModel() { return m_lat_model; };
	const LatticeDistanceModel& putLatticeDistanceModel() const { return m_lat_model; };
	const IPeakShiftFunc& putPeakShiftFunction() const { return *m_peak_shift_model; };
};

#endif

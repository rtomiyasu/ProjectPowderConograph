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
#ifndef _LatticeDistanceModel_hh_
#define _LatticeDistanceModel_hh_

// LatticeDistanceModel.hh

#include "../utility_data_structure/SymMat.hh"
#include "../utility_data_structure/VecDat3.hh"
#include "../zparam/ZParawError.hh"
#include "IMarquardtFmodel.hh"


// Interface class for LemarqMethod class.

class LatticeDistanceModel : public IMarquardtFmodel< VecDat3<Int4> >
{
private:
	enum{ ParaNum = 6 };
//	static bool m_param_view_flag;

	SymMat<ZParawError> m_param;	// s11, s22, s33, s23, s13, s12.
	Int4 m_num_indep;
	Mat_DP_constr m_cmat;

protected:
    Int4 setParam(Double*, ostringstream& outPutString);

	// Return the value(Ymod) of the model function at X.
	// If dYdX is not NULL, the entries of dYdP are the derivatives at X by the fitting parameters.
	bool putFittingCoef(const VecDat3<Int4>& X, Double& Ymod, Double* dYdP, ostringstream& outPutString) const;

	
public:
	LatticeDistanceModel();
	virtual ~LatticeDistanceModel();
	
	Int4 putParamNum() const;
	
	void putResult(SymMat<Double>& ans) const;
	
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
};

#endif

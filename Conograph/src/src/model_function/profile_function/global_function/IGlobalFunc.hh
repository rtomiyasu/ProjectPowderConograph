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
#ifndef _IGlobalFunc_hh_
#define _IGlobalFunc_hh_

// IGlobalFunc.hh
#include "../../../zparam/ZParawError.hh"
#include "../../../symmetric_operation/MillerIndex.hh"
#include "../../IMarquardtFmodel.hh"


// Interface class of global model functions for PhaseModelFunction class.
class IGlobalFunc : public IMarquardtFmodel<Double>
{

protected:
	// Return the value(Ymod) of the model function at X.
	// If dYdX is not NULL, also return the derivatives at X by dYdX.
	virtual bool putFittingCoef(const Double& X, Double& Ymod, Double* dYdX, ostringstream& outPutString) const = 0;

	
public:
    virtual ~IGlobalFunc(){};

	virtual MillerIndex3 putGlbAnisoHKL() const { return MillerIndex3(); };
	
	// Return the values and the standard errors of the parameters.
	virtual void putResult(vector<ZParawError>&) const = 0;

	// Return the number of parameters.( All the parameters are independently fit. )
    Int4 putNumberOfIndependentParam() const;

	// Returns the flags whether the parameter is linear or not.
//	virtual void putLinearParamFlag(bool* ans) const{ for(Int4 i=0; i<putNumberOfIndependentParam(); i++) ans[i] = false; };

	// Return the miller index that represents anisotropy.
	virtual void putDefaultFitFlag(vector<etype_ID>&) const;

	// Returns the flags on fix/fit.
	virtual void putFitFlag(vector<etype_ID>&) const = 0;

	// Set the constraints which indicates which parameters are fixed or independently fit.
	virtual void setConstraint(const constr_DP*) = 0;
	
	// Return constraints that indicates all the isotropic parameters are independently fit. 
	void putConstraint(constr_DP* const) const;	// This function always returns the constraint of isotropy. 

	// Return the values of the pararmeters that are fit to the local parameters.
	void setFittedParam(const Vec_DP&, const vector<Double>&, const vector<Double>&, 
							const vector<bool>& nxfit,
							const bool& output_view_flag,
							const Double&,  const Int4&, const vector<Double>&, Double& chisq_all);

};


#endif

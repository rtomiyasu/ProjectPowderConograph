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
#include"IGlobalFunc.hh"
#include"../../../levenberg_marquardt/LemarqMethod.hh"
#include"../../MarquardtFmodelBase.hh"


Int4 IGlobalFunc::putNumberOfIndependentParam() const
{
	Mat_DP_constr cmat( putParamNum() );
	putConstraint(&cmat[0]);
	return countNumberOfIndependentParam(cmat.begin(),cmat.end());
}

// Return the miller index that represents anisotropy.
void IGlobalFunc::putDefaultFitFlag(vector<etype_ID>& flag) const
{
	flag.clear();
	flag.resize( putParamNum(), _ZRietveldIDVary );
}

void IGlobalFunc::putConstraint(constr_DP* const tray) const
{
	vector<etype_ID> fitflag;
	putFitFlag( fitflag );
	
	for(Int4 k=0; k<putParamNum(); k++)
	{
		if( fitflag[k] == _ZRietveldIDFixed ) tray[k].ID = _ZRietveldIDFixed;
		else tray[k].ID = _ZRietveldIDVary;
		tray[k].constr.clear();	// No constraints.
	}
}


// If wt != 2, lclsterr is used as weights in least square method. 
// If wt = 2, lclsterr is not used,  all the weights equal 1. 
void IGlobalFunc::setFittedParam(const Vec_DP& dwidth, 
		const vector<Double>& lclparam, const vector<Double>& lclsterr,
		const vector<bool>& nxfit,
		const bool& output_view_flag, const Double& judge_conv,  const Int4& max_itnum,
		const vector<Double>& glbinit, Double& chisq_all)
{
	Int4 itnum;
	
	LemarqMethod marq;
	marq.setParam(output_view_flag, judge_conv, max_itnum, 0.001);
	
	MarquardtFmodelBase<Double> LMFunc(*this);

	LMFunc.addHistogram(dwidth, lclparam, lclsterr, nxfit);
	
	Vec_DP chisq;
	Int4 reason_terminate_LM;
	marq.execute(glbinit, LMFunc, itnum, chisq_all, chisq, reason_terminate_LM);
}

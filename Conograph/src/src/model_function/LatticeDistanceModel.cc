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
#include "LatticeDistanceModel.hh"
#include "../levenberg_marquardt/SVdcmp.hh"
#include "../levenberg_marquardt/LemarqMethod.hh"

LatticeDistanceModel::LatticeDistanceModel()
	: m_param(3, 0.0)
{
	m_cmat.resize(ParaNum);
	for(Int4 i=0; i<ParaNum; i++) m_cmat[i].ID = _ZRietveldIDFixed;
	m_num_indep = 0;
}

LatticeDistanceModel::~LatticeDistanceModel()
{
}

Int4 LatticeDistanceModel::putParamNum() const
{ 
	return ParaNum;
};

void LatticeDistanceModel::putResult(SymMat<Double>& ans) const
{
	for(Int4 i=0; i<3; i++)
	{
		for(Int4 j=0; j<=i; j++) ans(i,j) = m_param(i,j).value;
	}
}


Int4 LatticeDistanceModel::putNumberOfIndependentParam() const
{
	return m_num_indep;
}


void LatticeDistanceModel::putConstraint(constr_DP* const tray) const
{
	for(Int4 k=0; k<ParaNum; k++)
	{
		tray[k] = m_cmat[k];
	}
}



void LatticeDistanceModel::setConstraint(const constr_DP* tray)
{
	for(Int4 k=0; k<ParaNum; k++)
	{
		m_cmat[k] = tray[k];
	}
	m_num_indep = countNumberOfIndependentParam(m_cmat.begin(), m_cmat.end());
}


void LatticeDistanceModel::setCovariantMatrixAll(const SymMat<Double>& cov)
{
	m_param(0,0).error = sqrt_d( cov(0,0) );
	m_param(1,1).error = sqrt_d( cov(1,1) );
	m_param(2,2).error = sqrt_d( cov(2,2) );
	m_param(1,2).error = sqrt_d( cov(3,3) );
	m_param(0,2).error = sqrt_d( cov(4,4) );
	m_param(0,1).error = sqrt_d( cov(5,5) );
}


// Check if the symmetric matrix U is positive definite. 
Int4 checkSymMat3SemipositiveDefinite(const SymMat<ZParawError>& U)
{
	NRMat<Double> U_mat(3,3);

	// Set U_mat.
	for(Int4 i=0; i<3; i++) U_mat[i][i] = U(i,i).value;
    U_mat[0][1] = U(0,1).value;
    U_mat[1][0] = U(0,1).value;
    U_mat[0][2] = U(0,2).value;
    U_mat[2][0] = U(0,2).value;
    U_mat[1][2] = U(1,2).value;
    U_mat[2][1] = U(1,2).value;
			
	// SV decomposition.
    NRVec<Double> d(3);
	SVdcmp(U_mat, d, 0, 3);

	Int4 ans = 1;
	for(Int4 i=0; i<3; i++)
	{
		if(d[i] < 0.0)
		{
			d[i] = 0.0;
			ans = 2;
		}
	}

	return ans;
}


Int4 LatticeDistanceModel::setParam(Double* param, ostringstream& outPutString)
{
	m_param(0,0) = param[0];
	m_param(1,1) = param[1];
	m_param(2,2) = param[2];
	m_param(1,2) = param[3];
	m_param(0,2) = param[4];
	m_param(0,1) = param[5];
	return checkSymMat3SemipositiveDefinite(m_param);
}


bool LatticeDistanceModel::putFittingCoef(const VecDat3<Int4>& hkl, Double& lclparam, Double* dlclparam_dother, ostringstream& outPutString) const
{
	lclparam = 0.0;
	const Double hkl00 = hkl[0]*hkl[0];
	const Double hkl11 = hkl[1]*hkl[1];
	const Double hkl22 = hkl[2]*hkl[2];
	const Double hkl01 = hkl[0]*hkl[1]*2;
	const Double hkl02 = hkl[0]*hkl[2]*2;
	const Double hkl12 = hkl[1]*hkl[2]*2;

	lclparam = m_param(0,0).value*hkl00;
	lclparam += m_param(1,1).value*hkl11;
	lclparam += m_param(2,2).value*hkl22;
	lclparam += m_param(0,1).value*hkl01;
	lclparam += m_param(0,2).value*hkl02;
	lclparam += m_param(1,2).value*hkl12;
	
	if( dlclparam_dother != NULL )
	{
		Int4 j;
		for(j=0; j<m_num_indep; j++) dlclparam_dother[j] = 0.0;

		j=0;
		if( m_cmat[0].ID == _ZRietveldIDVary ) dlclparam_dother[j++]+=hkl00;
		else if( m_cmat[0].ID == _ZRietveldIDDepend )
		{
			for(Vec_DP_save::const_iterator it=m_cmat[0].constr.begin(); it<m_cmat[0].constr.end(); it++)
			{
				dlclparam_dother[it->index] += it->element * hkl00;
			}
		}

		if( m_cmat[1].ID == _ZRietveldIDVary ) dlclparam_dother[j++]+=hkl11;
		else if( m_cmat[1].ID == _ZRietveldIDDepend )
		{
			for(Vec_DP_save::const_iterator it=m_cmat[1].constr.begin(); it<m_cmat[1].constr.end(); it++)
			{
				dlclparam_dother[it->index] += it->element * hkl11;
			}
		}

		if( m_cmat[2].ID == _ZRietveldIDVary ) dlclparam_dother[j++]+=hkl22;
		else if( m_cmat[2].ID == _ZRietveldIDDepend )
		{
			for(Vec_DP_save::const_iterator it=m_cmat[2].constr.begin(); it<m_cmat[2].constr.end(); it++)
			{
				dlclparam_dother[it->index] += it->element * hkl22;
			}
		}

		if( m_cmat[3].ID == _ZRietveldIDVary ) dlclparam_dother[j++]+=hkl12;
		else if( m_cmat[3].ID == _ZRietveldIDDepend )
		{
			for(Vec_DP_save::const_iterator it=m_cmat[3].constr.begin(); it<m_cmat[3].constr.end(); it++)
			{
				dlclparam_dother[it->index] += it->element * hkl12;
			}
		}

		if( m_cmat[4].ID == _ZRietveldIDVary ) dlclparam_dother[j++]+=hkl02;
		else if( m_cmat[4].ID == _ZRietveldIDDepend )
		{
			for(Vec_DP_save::const_iterator it=m_cmat[4].constr.begin(); it<m_cmat[4].constr.end(); it++)
			{
				dlclparam_dother[it->index] += it->element * hkl02;
			}
		}

		if( m_cmat[5].ID == _ZRietveldIDVary ) dlclparam_dother[j++]+=hkl01;
		else if( m_cmat[5].ID == _ZRietveldIDDepend )
		{
			for(Vec_DP_save::const_iterator it=m_cmat[5].constr.begin(); it<m_cmat[5].constr.end(); it++)
			{
				dlclparam_dother[it->index] += it->element * hkl01;
			}
		}
	}

	return true;
}


// If wt != 2, lclsterr is used as weights in least square method. 
// If wt = 2, lclsterr is not used,  all the weights equal 1. 
pair<bool, ZErrorMessage> LatticeDistanceModel::setFittedParam(const vector< VecDat3<Int4> >& hkl,
		const vector<Double>& lclparam, const vector<Double>& lclvar,
		const vector<bool>& nxfit,
		const bool& output_view_flag, const Double& judge_conv,  const Int4& max_itnum,
		const vector<Double>& glbinit, Double& chisq_all)
{
	Int4 itnum;
	
	LemarqMethod marq;
	marq.setParam(output_view_flag, judge_conv, max_itnum, 0.001);
	
	MarquardtFmodelBase< VecDat3<Int4> > LMFunc(*this);
	ZErrorMessage zerr = LMFunc.addHistogram(hkl, lclparam, lclvar, nxfit);
	if( zerr.putErrorType() != ZErrorNoError ) return pair<bool, ZErrorMessage>(false, zerr);
	
	Vec_DP chisq;
	Int4 reason_terminate_LM;
	return marq.execute< VecDat3<Int4> >(glbinit, LMFunc, itnum, chisq_all, chisq, reason_terminate_LM);
}

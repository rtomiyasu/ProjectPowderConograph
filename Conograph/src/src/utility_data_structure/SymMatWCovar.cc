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
#include "SymMatWCovar.hh"

SymMatWCovar::SymMatWCovar(const Int4& isize)
: ISIZE(isize), m_S(isize), m_S_covar(isize*(isize+1)/2)
{
}


SymMatWCovar::SymMatWCovar(const SymMat<Double>& S, const SymMat<Double>& S_covar)
	: ISIZE(S.size()), m_S(S), m_S_covar(S_covar)
{
	assert( ISIZE*(ISIZE+1)/2 == S_covar.size() ); 
}


SymMatWCovar::SymMatWCovar(const SymMatWCovar& rhs)
	: ISIZE(rhs.ISIZE), m_S(rhs.m_S), m_S_covar(rhs.m_S_covar)
{
}


SymMatWCovar::~SymMatWCovar()
{
}


SymMatWCovar& SymMatWCovar::operator=(const SymMatWCovar &rhs)
{
	if (this != &rhs)
	{
		assert( ISIZE == rhs.ISIZE );
		m_S = rhs.m_S;
		m_S_covar = rhs.m_S_covar;
	}
	return *this;
}


Double SymMatWCovar::Determinant3(Double& var) const
{
	assert( ISIZE == 3 );
	
	Double det_coef[6];	// 6 = 3*(3+1)/2
	
	det_coef[0] = m_S(1,1)*m_S(2,2)-m_S(1,2)*m_S(1,2);	// = d(detS)/d(m_S[0])
	det_coef[1] = - (m_S(0,1)*m_S(2,2)-m_S(0,2)*m_S(1,2));	// = 0.5 * d(detS)/d(m_S[1])
	det_coef[2] = m_S(0,1)*m_S(1,2)-m_S(0,2)*m_S(1,1);	// = 0.5 * d(detS)/d(m_S[2])
	det_coef[3] = m_S(0,0)*m_S(2,2)-m_S(0,2)*m_S(0,2);	// = d(detS)/d(m_S[3])
	det_coef[4] = - (m_S(0,0)*m_S(1,2)-m_S(0,1)*m_S(0,2));	// = 0.5 * d(detS)/d(m_S[4])
	det_coef[5] = m_S(0,0)*m_S(1,1)-m_S(0,1)*m_S(0,1);	// = d(detS)/d(m_S[5])

	const Double det_rhs = m_S(0,0)*det_coef[0] + m_S(0,1)*det_coef[1] + m_S(0,2)*det_coef[2];

	det_coef[1] *= 2.0;
	det_coef[2] *= 2.0;
	det_coef[4] *= 2.0;

	var = 0.0;
	for(Int4 k=0; k<6; k++)
	{
		for(Int4 k2=0; k2<k; k2++) var += 2.0 * m_S_covar(k, k2) * (det_coef[k] * det_coef[k2]);
		var += m_S_covar(k, k) * det_coef[k] * det_coef[k];
	}
	
	return det_rhs;
}



SymMatWCovar SymMatWCovar::Inverse3() const 
{
	assert( this->size() == 3 );
	
	const Double det12 = m_S(0,0)*m_S(1,1)-m_S(0,1)*m_S(1,0);	// = d(detS)/dC*
	const Double det13 = m_S(0,0)*m_S(2,2)-m_S(0,2)*m_S(2,0);	// = d(detS)/dB*
	const Double det23 = m_S(1,1)*m_S(2,2)-m_S(1,2)*m_S(2,1);	// = d(detS)/dA*
	const Double det23_13 = m_S(1,0)*m_S(2,2)-m_S(1,2)*m_S(2,0);	// = - 0.5 * d(detS)/dF*
	const Double det23_12 = m_S(1,0)*m_S(2,1)-m_S(1,1)*m_S(2,0);	// = 0.5 * d(detS)/dE*
	const Double det13_12 = m_S(0,0)*m_S(2,1)-m_S(0,1)*m_S(2,0);	// = - 0.5 * d(detS)/dD*

	const Double det = m_S(0,0)*det23 - m_S(0,1)*det23_13 + m_S(0,2)*det23_12;
	const Double inv_det = 1.0/det;
	
	SymMatWCovar ans(3);
	ans.m_S(0,0) = det23;	// ans(0,0)
	ans.m_S(1,1) = det13;	// ans(1,1)
	ans.m_S(2,2) = det12;	// ans(2,2)
	ans.m_S(0,1) = -det23_13;	// ans(0,1)
	ans.m_S(0,2) = det23_12;	// ans(0,2)
	ans.m_S(1,2) = -det13_12;	// ans(1,2)
	ans.m_S *= inv_det;

	// Calculate inv_covar.
	NRMat<Double> cmat(6,6);
	// det23 = m_S(1,1)*m_S(2,2)-m_S(1,2)*m_S(2,1).
	cmat[0][0] = - det23*det23;	// dA / dA*
	cmat[0][1] = + 2.0*det23*det23_13;	// dA / dF*
	cmat[0][2] = - 2.0*det23*det23_12;	// dA / dE*
	cmat[0][3] = - det23*det13 + m_S(2,2)*det;	// dA / dB*
	cmat[0][4] = + 2.0*det23*det13_12 - 2.0*m_S(1,2)*det;	// dA / dD*
	cmat[0][5] = - det23*det12 + m_S(1,1)*det;	// dA / dC*

	// -det23_13 = -m_S(2,2)*m_S(1,0)+m_S(1,2)*m_S(2,0).
	cmat[1][0] = + det23_13*det23 ;	// dF / dA*
	cmat[1][1] = - 2.0*det23_13*det23_13 - m_S(2,2)*det;	// dF / dF*
	cmat[1][2] = + 2.0*det23_13*det23_12 + m_S(1,2)*det;	// dF / dE*
	cmat[1][3] = + det23_13*det13;	// dF / dB*
	cmat[1][4] = - 2.0*det23_13*det13_12 + m_S(2,0)*det;	// dF / dD*
	cmat[1][5] = + det23_13*det12 - m_S(1,0)*det;	// dF / dC*

	// det23_12 = m_S(1,0)*m_S(2,1)-m_S(2,0)*m_S(1,1).
	cmat[2][0] = - det23_12*det23;	// dE / dA*
	cmat[2][1] = + 2.0*det23_12*det23_13 + m_S(2,1)*det;	// dE / dF*
	cmat[2][2] = - 2.0*det23_12*det23_12 - m_S(1,1)*det;	// dE / dE*
	cmat[2][3] = - det23_12*det13 - m_S(2,0)*det;	// dE / dB*
	cmat[2][4] = + 2.0*det23_12*det13_12 + m_S(1,0)*det;	// dE / dD*
	cmat[2][5] = - det23_12*det12;	// dE / dC*

	// det13 = m_S(0,0)*m_S(2,2)-m_S(0,2)*m_S(2,0).
	cmat[3][0] = - det13*det23 + m_S(2,2)*det;	// dB / dA*
	cmat[3][1] = + 2.0*det13*det23_13;	// dB / dF*
	cmat[3][2] = - 2.0*det13*det23_12 - 2.0*m_S(0,2)*det;	// dB / dE*
	cmat[3][3] = - det13*det13;	// dB / dB*
	cmat[3][4] = + 2.0*det13*det13_12;	// dB / dD*
	cmat[3][5] = - det13*det12 + m_S(0,0)*det;	// dB / dC*

	//	-det13_12 = -m_S(0,0)*m_S(2,1)+m_S(0,1)*m_S(2,0).
	cmat[4][0] = + det13_12*det23 - m_S(2,1)*det;	// dD / dA*
	cmat[4][1] = - 2.0*det13_12*det23_13 + m_S(2,0)*det;	// dD / dF*
	cmat[4][2] = + 2.0*det13_12*det23_12 + m_S(0,1)*det;	// dD / dE*
	cmat[4][3] = + det13_12*det13;	// dD / dB*
	cmat[4][4] = - 2.0*det13_12*det13_12 - m_S(0,0)*det;	// dD / dD*
	cmat[4][5] = + det13_12*det12;	// dD / dC*

	// det12 = m_S(0,0)*m_S(1,1)-m_S(0,1)*m_S(1,0).
	cmat[5][0] = - det12*det23 + m_S(1,1)*det;	// dC / dA*
	cmat[5][1] = + 2.0*det12*det23_13 - 2.0*m_S(0,1)*det;	// dC / dF*
	cmat[5][2] = - 2.0*det12*det23_12;	// dC / dE*
	cmat[5][3] = - det12*det13 + m_S(0,0)*det;	// dC / dB*
	cmat[5][4] = + 2.0*det12*det13_12;	// dC / dD*
	cmat[5][5] = - det12*det12;	// dC / dC*

	// inv_covar = cmat * m_S_covar * transpose(cmat).
	Double inv_detsq = inv_det*inv_det;
	inv_detsq *= inv_detsq;

	ans.m_S_covar = transform_sym_matrix( cmat, m_S_covar );

	ans.m_S_covar *= inv_detsq;

	return ans;
}

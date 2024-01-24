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
#ifndef SymMatWCovar_HH_
#define SymMatWCovar_HH_

#include "../RietveldAnalysisTypes.hh"
#include "../utility_data_structure/SymMat.hh"
#include "../utility_data_structure/nrutil_nr.hh"
#include "../utility_func/transform_sym_matrix.hh"

class SymMatWCovar
{
	friend inline SymMatWCovar transform_sym_matrix(const NRMat<Int4>& lhs, const SymMatWCovar& rhs);

private:
	const Int4 ISIZE;
	SymMat<Double> m_S;
	SymMat<Double> m_S_covar;
	
	SymMatWCovar(const SymMat<Double>& S, const SymMat<Double>& S_covar);

public:
	SymMatWCovar(const Int4& isize);
	SymMatWCovar(const SymMatWCovar& rhs);
	~SymMatWCovar();

	SymMatWCovar& operator=(const SymMatWCovar &rhs);	//assignment

	inline const Int4& size() const { return ISIZE; };

	inline Double& operator()(const Int4& i, const Int4& j) { return m_S(i,j); };
	inline const Double& operator()(const Int4& i, const Int4& j) const { return m_S(i,j); };
	inline SymMatWCovar& operator*=(const Double& t){ m_S *= t; m_S_covar *= t*t; return *this; };
	
	inline const Double& Var(const Int4& i, const Int4& j) const { const Int4 index = put_index(ISIZE,i,j); return m_S_covar(index, index); };
	inline Double& Covar(const Int4& i, const Int4& j, const Int4& k, const Int4& l) { return m_S_covar(put_index(ISIZE,i,j), put_index(ISIZE,k,l)); };
	inline const Double& Covar(const Int4& i, const Int4& j, const Int4& k, const Int4& l) const { return m_S_covar(put_index(ISIZE,i,j), put_index(ISIZE,k,l)); };
	
	inline const SymMat<Double>& ValMat() const { return m_S; };
	inline const SymMat<Double>& CovMat() const { return m_S_covar; };

	Double Determinant3(Double& var) const;
	SymMatWCovar Inverse3() const;
};


inline SymMatWCovar transform_sym_matrix(const NRMat<Int4>& lhs, const SymMatWCovar& rhs)
{ 
	return SymMatWCovar( transform_sym_matrix(lhs, rhs.m_S), transform_sym_covar(lhs, rhs.m_S_covar) );
};


#endif /*SymMatWCovar_HH_*/

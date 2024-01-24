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
#ifndef _Node3_HH_
#define _Node3_HH_

#include <assert.h>
#include "../RietveldAnalysisTypes.hh"
#include "../zerror_type/error_out.hh"
#include "../utility_func/chToDouble.hh"

// (0,1) -> 0, (0,2) -> 1, (1,2) -> 2,
inline Int4 put_index_ij(const Int4& i, const Int4& j)
{
	assert( 0 <= i && 0 <= j && i != j && i < 3 && j < 3 );
	return i + j - 1;
}

class Node3
{
private:
	enum{ NUM_K = 3, NUM_SUM_K = 3 };
	VCData m_K[NUM_K];
	VCData m_sum_K[NUM_SUM_K];

	Double m_detS;
	Double m_detS_var;

	inline VCData& K(const Int4& i);
	inline VCData& SumK(const Int4& i, const Int4& j);

public:
	inline const VCData& K(const Int4& i) const;
	inline const VCData& SumK(const Int4& i, const Int4& j) const;

	inline Int4 iK(const Int4& i) const;
	inline Int4 iSumK(const Int4& i, const Int4& j) const;

	bool setBud0(const Double& cv2,
					const VCData& K0,
					const VCData& K1,
					const VCData& K2,
					const VCData& K3,
					const VCData& K01,
					const VCData& K02,
					SymMat<VCData>& S);

	Node3();
	Node3(const Node3&);
	Node3(const Double& arg){ m_detS = arg; };
	~Node3();
	
	Node3& operator=(const Node3&);
	inline bool operator<(const Node3&) const;
	inline bool operator==(const Node3&) const;

	inline const Double& putDeterminantS() const { return m_detS; }; 
};


inline VCData& Node3::K(const Int4& i)
{
	assert( 0 <= i && i < NUM_K );
	return m_K[i]; 
};

inline const VCData& Node3::K(const Int4& i) const
{
	assert( 0 <= i && i < NUM_K );
	return m_K[i];
};


inline Int4 Node3::iK(const Int4& i) const
{
	assert( 0 <= i && i < NUM_K );
	return m_K[i].putIndex();
};



inline VCData& Node3::SumK(const Int4& i, const Int4& j)
{
	return m_sum_K[ put_index_ij(i, j) ];
};


inline const VCData& Node3::SumK(const Int4& i, const Int4& j) const
{
	return m_sum_K[ put_index_ij(i, j) ];
};


// (0,1) -> 0, (0,2) -> 1, (1,2) -> 2, (0,3) -> 2, (1,3) -> 1, (2,3) -> 0.
inline Int4 Node3::iSumK(const Int4& i, const Int4& j) const
{
	return m_sum_K[ put_index_ij(i, j) ].putIndex();
};


inline bool Node3::operator<(const Node3& rhs) const
{
	return m_detS < rhs.m_detS;
}


inline Double Determinant3(const SymMat<VCData>& rhs, Double& var)
{
	return chToDoubleWCovar(rhs).Determinant3(var);
}


inline SymMat<VCData> calculateS3(const Node3& nodex)
{
	SymMat<VCData> ans(3);
	
	ans(0,0) = nodex.K(0);
	ans(0,1) = ( nodex.SumK(0,1) - nodex.K(0) - nodex.K(1) ) / 2;
	ans(0,2) = ( nodex.SumK(0,2) - nodex.K(0) - nodex.K(2) ) / 2;
	ans(1,1) = nodex.K(1);
	ans(1,2) = ( nodex.SumK(1,2) - nodex.K(1) - nodex.K(2) ) / 2;
	ans(2,2) = nodex.K(2);
	return ans;
}

#endif

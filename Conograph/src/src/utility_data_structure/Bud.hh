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
#ifndef _Bud_hh_
#define _Bud_hh_

#include <set>
#include "../RietveldAnalysisTypes.hh"
#include "SymMat.hh"
#include "VCData.hh"

// Bud.hh

class Bud
{
	friend inline bool operator==(const Bud& lhs, const Bud& rhs);
	friend inline bool operator<(const Bud& lhs, const Bud& rhs);
	friend inline Bud putBud124(const Int4& K1, const Int4& K2, Int4 K4=-1);
	friend inline Bud putBud123(const Int4& K1, const Int4& K2, const Int4& K3);

private:
	Int4 m_K[4];	// K[0] <= K[1], Q[2] <= Q[3].

public:

	Bud();
//	Bud(const Int4& K1, const Int4& K2, const Int4& K3, const Int4& K4, bool swap_flag=true);
	~Bud();

	void setIndex(const Int4& K1, const Int4& K2, const Int4& K3, const Int4& K4);

	inline const Int4& iK1() const { return m_K[0]; };
	inline const Int4& iK2() const { return m_K[1]; };
	inline const Int4& iK3() const { return m_K[2]; };
	inline const Int4& iK4() const { return m_K[3]; };

	inline VCData Q1() const;
	inline VCData Q2() const;
	inline VCData Q3() const;
	inline VCData Q4() const;

	inline bool IsSuperBasisObtuse() const;

	inline Double inner_product_312() const;
	inline Double cross_product_312() const;
	inline Double cross_product_412() const;
};

inline VCData Bud::Q1() const
{
	if( m_K[0] >= 0 ) return VCData(m_K[0], 1);
	else return ( VCData(m_K[2], 1) + VCData(m_K[3], 1) ) / 2 - VCData(m_K[1], 1);
};


inline VCData Bud::Q2() const
{
	if( m_K[1] >= 0 ) return VCData(m_K[1], 1);
	else return ( VCData(m_K[2], 1) + VCData(m_K[3], 1) ) / 2 - VCData(m_K[0], 1);
};

inline VCData Bud::Q3() const
{
	if( m_K[2] >= 0 ) return VCData(m_K[2], 1);
	else return ( VCData(m_K[0], 2) + VCData(m_K[1], 2) ) * 2 - VCData(m_K[3], 1);
};

inline VCData Bud::Q4() const
{
	if( m_K[3] >= 0 ) return VCData(m_K[3], 1);
	else return ( VCData(m_K[0], 2) + VCData(m_K[1], 2) ) * 2 - VCData(m_K[2], 1);
};

inline bool Bud::IsSuperBasisObtuse() const
{
	Double Q1, Q2, Q3;
	Q1 = this->Q1().Value();
	Q2 = this->Q2().Value();
	Q3 = this->Q3().Value();

	if( Q2 + Q3 <= Q1 ) return false; 
	if( Q1 + Q3 <= Q2 ) return false; 
	if( Q1 + Q2 <= Q3 ) return false; 
 	else return true;
}

inline Double Bud::inner_product_312() const
{
	Double Q1, Q2, Q3;
	Q1 = this->Q1().Value();
	Q2 = this->Q2().Value();
	Q3 = this->Q3().Value();

	return ( Q3 - Q1 - Q2 )*0.5;
}

inline Double Bud::cross_product_312() const
{
	Double Q1, Q2, Q3;
	Q1 = this->Q1().Value();
	Q2 = this->Q2().Value();
	Q3 = this->Q3().Value();
	
	Double inn = ( Q3 - Q1 - Q2 )*0.5;

	return Q1*Q2 - inn*inn;
}

inline Double Bud::cross_product_412() const
{
	Double Q1, Q2, Q4;
	Q1 = this->Q1().Value();
	Q2 = this->Q2().Value();
	Q4 = this->Q4().Value();

	Double inn = ( Q4 - Q1 - Q2 )*0.5;

	return Q1*Q2 - inn*inn;
}


inline bool operator==(const Bud& lhs, const Bud& rhs)
{
	if( lhs.m_K[0] != rhs.m_K[0] ) return false;
	if( lhs.m_K[1] != rhs.m_K[1] ) return false;
	if( lhs.m_K[2] != rhs.m_K[2] ) return false;
	if( lhs.m_K[3] != rhs.m_K[3] ) return false;
	return true;
}

inline bool operator<(const Bud& lhs, const Bud& rhs)
{
	if( lhs.m_K[0] < rhs.m_K[0] ) return true;
	if( lhs.m_K[0] > rhs.m_K[0] ) return false;

	if( lhs.m_K[1] < rhs.m_K[1] ) return true;
	if( lhs.m_K[1] > rhs.m_K[1] ) return false;

	if( lhs.m_K[3] < rhs.m_K[3] ) return true;
	if( lhs.m_K[3] > rhs.m_K[3] ) return false;

	if( lhs.m_K[2] < rhs.m_K[2] ) return true;
	return false;
}


inline Bud putBud124(const Int4& K1, const Int4& K2, Int4 K4)
{
	Bud budex;
	if( K1 <= K2 )
	{
		budex.m_K[0] = K1;
		budex.m_K[1] = K2;
	}
	else
	{
		budex.m_K[0] = K2;
		budex.m_K[1] = K1;
	}
	budex.m_K[2] = -1;
	budex.m_K[3] = K4;
	return budex;
}


inline Bud putBud123(const Int4& K1, const Int4& K2, const Int4& K3)
{
	Bud budex;
	if( K1 <= K2 )
	{
		budex.m_K[0] = K1;
		budex.m_K[1] = K2;
	}
	else
	{
		budex.m_K[0] = K2;
		budex.m_K[1] = K1;
	}
	budex.m_K[2] = K3;
	budex.m_K[3] = -1;
	return budex;
}

template<class Iterator>
inline pair< Iterator, Iterator> equal_range_bud12(
		const Iterator& it_begin,
		const Iterator& it_end, 
		const Int4& K1, const Int4& K2)
{
	pair<Iterator, Iterator> it_pair;
	const Bud budex = putBud124(K1, K2);
	it_pair.first = lower_bound(it_begin, it_end, budex);
	it_pair.second = lower_bound(it_pair.first, it_end, putBud124(budex.iK1(), budex.iK2()+1));

	return it_pair;
}

template<class Iterator>
inline pair< Iterator, Iterator> equal_range_bud124(
		const Iterator& it_begin,
		const Iterator& it_end, 
		const Int4& K1, const Int4& K2, const Int4& K4)
{
	pair<Iterator, Iterator> it_pair;
	const Bud budex = putBud124(K1, K2, K4);
	it_pair.first = lower_bound(it_begin, it_end, budex);
	it_pair.second = lower_bound(it_pair.first, it_end, putBud124(budex.iK1(), budex.iK2(), budex.iK4()+1));

	return it_pair;
}


#endif

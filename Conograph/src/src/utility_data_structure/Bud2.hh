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
#ifndef _Bud2_hh_
#define _Bud2_hh_

#include <set>
#include "../RietveldAnalysisTypes.hh"
#include "SymMat.hh"
#include "VCData.hh"

// Bud2.hh

class Bud2
{
	friend inline bool operator==(const Bud2& lhs, const Bud2& rhs);
	friend inline bool operator<(const Bud2& lhs, const Bud2& rhs);

private:
	Int4 m_K[4];

public:

	Bud2();
	~Bud2();

	void setIndex(const Int4& K1, const Int4& K2, const Int4& K3, const Int4& K4);

	inline const Int4& iK1() const { return m_K[0]; };
	inline const Int4& iK2() const { return m_K[1]; };
	inline const Int4& iK12() const { return m_K[2]; };
	inline const Int4& iK1_2() const { return m_K[3]; };

	inline VCData Q1() const;
	inline VCData Q2() const;
	inline VCData Q12() const;
//	inline VCData Q1_2() const;

	static Bud2 putBud2(const Int4& K1, const Int4 K2=-1, const Int4 K3=-1, const Int4 K4=-1)
	{
		Bud2 nodex;
		nodex.m_K[0] = K1;
		nodex.m_K[1] = K2;
		nodex.m_K[2] = K3;
		nodex.m_K[3] = K4;

		return nodex;
	}
};

inline VCData Bud2::Q1() const
{
	if( m_K[0] >= 0 ) return VCData(m_K[0], 1);
	else return ( VCData(m_K[2], 1) + VCData(m_K[3], 1) ) / 2 - VCData(m_K[1], 1);
};


inline VCData Bud2::Q2() const
{
	if( m_K[1] >= 0 ) return VCData(m_K[1], 1);
	else return ( VCData(m_K[2], 1) + VCData(m_K[3], 1) ) / 2 - VCData(m_K[0], 1);
};

inline VCData Bud2::Q12() const
{
	if( m_K[2] >= 0 ) return VCData(m_K[2], 1);
	else return ( VCData(m_K[0], 2) + VCData(m_K[1], 2) ) * 2 - VCData(m_K[3], 1);
};

//inline VCData Bud2::Q1_2() const
//{
//	if( m_K[3] >= 0 ) return VCData(m_K[3], 1);
//	else return ( VCData(m_K[0], 2) + VCData(m_K[1], 2) ) * 2 - VCData(m_K[2], 1);
//};


inline bool operator==(const Bud2& lhs, const Bud2& rhs)
{
	if( lhs.m_K[0] != rhs.m_K[0] ) return false;
	if( lhs.m_K[1] != rhs.m_K[1] ) return false;
	if( lhs.m_K[2] != rhs.m_K[2] ) return false;
	if( lhs.m_K[0] >= 0 && lhs.m_K[1] >= 0 && lhs.m_K[2] >= 0 ) return true;
	if( lhs.m_K[3] != rhs.m_K[3] ) return false;
	return true;
}



inline bool operator<(const Bud2& lhs, const Bud2& rhs)
{
	if( lhs.m_K[0] < rhs.m_K[0] ) return true;
	if( lhs.m_K[0] > rhs.m_K[0] ) return false;

	if( lhs.m_K[1] < rhs.m_K[1] ) return true;
	if( lhs.m_K[1] > rhs.m_K[1] ) return false;

	if( lhs.m_K[2] < rhs.m_K[2] ) return true;
	if( lhs.m_K[2] > rhs.m_K[2] ) return false;
	if( lhs.m_K[0] >= 0 && lhs.m_K[1] >= 0 && lhs.m_K[2] >= 0 ) return false;

	if( lhs.m_K[3] < rhs.m_K[3] ) return true;
	return false;
}

#endif

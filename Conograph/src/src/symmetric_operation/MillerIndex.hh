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
#ifndef MILLERINDEX_HH_
#define MILLERINDEX_HH_

#include<complex>
#include<cmath>
#include"../RietveldAnalysisTypes.hh"
#include"../utility_data_structure/VecDat3.hh"
#include"../utility_data_structure/SymMat.hh"
#include"SymmetricOperation.hh"
#include"../utility_func/zmath.hh"

// class of hkl or hkli(h+k+i=0)
// In the rhombohedral case, m_hkl[3] is just a dummy. 
// In the hexagonal case(i.e. symmetric action given by symmetric operations belong to D6h), m_hkl[3] always equals to -m_hkl[0]-m_hkl[1].
class MillerIndex3
{
	friend inline MillerIndex3 operator*(const SymmetricOperation&, const MillerIndex3&);
	
private:
	enum{ISIZE=3};
	
	Int4 m_hkl[ISIZE+1]; // The last index used for calculation.( Supppose that m_hkl[3] = -m_hkl[1] - m_hkl[2] all the time.)
	inline Int4& operator[](const Int4&);
	
public:
	MillerIndex3();
	MillerIndex3(const MillerIndex3&);
	MillerIndex3(const VecDat3<Int4>&);
	MillerIndex3(const Int4&, const Int4&, const Int4&);
	~MillerIndex3();

	MillerIndex3& operator=(const MillerIndex3&);
	inline const Int4& operator[](const Int4&) const;
};

inline Int4& MillerIndex3::operator[](const Int4& k)
{
	return m_hkl[k];
}

inline const Int4& MillerIndex3::operator[](const Int4& k) const
{
	return m_hkl[k];
}

inline bool operator==(const MillerIndex3& lhs, const MillerIndex3& rhs)
{
	if(lhs[0] != rhs[0]) return false;
	if(lhs[1] != rhs[1]) return false;
	if(lhs[2] != rhs[2]) return false;
	return true;
}

//inline bool operator!=(const MillerIndex3& lhs, const MillerIndex3& rhs)
//{
//	return !(lhs==rhs);
//}

inline MillerIndex3 operator*(const SymmetricOperation& so, const MillerIndex3& hkl)
{
	MillerIndex3 hkl2;
	for(Int4 i=0; i<4; i++) hkl2[so.Perm[i]]=so.Sign[i]*hkl[i];
	return hkl2;
}

#endif /*MILLERINDEX_HH_*/


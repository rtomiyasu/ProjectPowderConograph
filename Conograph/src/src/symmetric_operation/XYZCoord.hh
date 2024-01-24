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
#ifndef XYZCoord_HH_
#define XYZCoord_HH_

#include "../RietveldAnalysisTypes.hh"
#include "StringS1.hh"
#include "SymmetricOperation.hh"

template <class T>
class XYZCoord
{
private:
	enum{ISIZE=3};     // The last index used for calculation.
	T m_xyz[ISIZE+1];	// m_xyz[3] equals alway1 0.
	inline T& operator[](const Int4&);

public:
	XYZCoord();
	XYZCoord(const XYZCoord<T>&);
	XYZCoord(const T&, const T&, const T&);
	~XYZCoord();

	XYZCoord<T>& operator=(const XYZCoord<T>&);
//	inline VecDat3<T> toVecDat3() const;

	template <class S>
	inline XYZCoord<T>& operator+=(const XYZCoord<S>&);
	template <class S>
	inline XYZCoord<T> operator+(const XYZCoord<S>&) const;
	
	inline XYZCoord<T> operator*(const SymmetricOperation&) const;

	inline const T& operator[](const Int4&) const;
};

template <class T>
XYZCoord<T>::XYZCoord()
{
	m_xyz[0]=0;
	m_xyz[1]=0;
	m_xyz[2]=0;
	m_xyz[3]=0;
}

template <class T>
XYZCoord<T>::XYZCoord(const XYZCoord<T>& rhs)
{
	for(Int4 k=0; k<ISIZE+1; k++) m_xyz[k] = rhs.m_xyz[k];
}

template <class T>
XYZCoord<T>::XYZCoord(const T& X, const T& Y, const T& Z)
{
	m_xyz[0]=X;
	m_xyz[1]=Y;
	m_xyz[2]=Z;
	m_xyz[3]=0;
}

template <class T>
XYZCoord<T>::~XYZCoord()
{
}

template <class T>
XYZCoord<T>& XYZCoord<T>::operator=(const XYZCoord<T>& rhs)
{
	if( this != &rhs)
		for(Int4 k=0; k<ISIZE+1; k++) m_xyz[k] = rhs.m_xyz[k];
	return *this;
}

/*
template <class T>
inline VecDat3<T> XYZCoord<T>::toVecDat3() const
{
	VecDat3<T> ans;
	for(Int4 k=0; k<ISIZE; k++) ans[k] = m_xyz[k];
	return ans;
}
*/

template <class T>
inline T& XYZCoord<T>::operator[](const int& k)
{
	return m_xyz[k];
}

template <class T>
inline const T& XYZCoord<T>::operator[](const int& k) const
{
	return m_xyz[k];
}

template <class T> template<class S>
inline XYZCoord<T>& XYZCoord<T>::operator+=(const XYZCoord<S>& rhs)
{
	for(int k=0; k<ISIZE+1; k++) m_xyz[k]+=rhs[k];
	return *this;
}

template <class T> template<class S>
inline XYZCoord<T> XYZCoord<T>::operator+(const XYZCoord<S>& rhs) const
{
	XYZCoord<T> ans(*this);
	for(int k=0; k<ISIZE+1; k++) ans[k]+=rhs[k];
	return ans;
}

template <class T>
inline XYZCoord<T> XYZCoord<T>::operator*(const SymmetricOperation& so) const
{
	XYZCoord<T> xyz2;
	for(Int4 i=0; i<4; i++)
	{
		if(so.Sign[i]>0) xyz2[i] = m_xyz[so.Perm[i]];
		else xyz2[i] = -m_xyz[so.Perm[i]];
	}
	xyz2[0]-=xyz2[3];
	xyz2[1]-=xyz2[3];
	xyz2[3]=0;
	
	return xyz2;
}

template <class T>
inline bool operator<(const XYZCoord<T>& lhs, const XYZCoord<T>& rhs)
{
	for(int k=0; k<3; k++)
	{
		if( lhs[k] < rhs[k] ) return true;
		if( rhs[k] < lhs[k] ) return false;
	}
	return false;
}

template <class T>
inline bool operator==(const XYZCoord<T>& lhs, const XYZCoord<T>& rhs)
{
	for(int k=0; k<3; k++)
		if( !(lhs[k] == rhs[k]) ) return false;
	return true;
}

// For rhs = (ax+by+cz, dx+ey+fz, gx+hy+iz), coef_mat is the square matrix of size 3 as follows on output.
//	a b c
//  d e f
//  g h i
inline void putMatrixForm(const XYZCoord<StringS1>& rhs, NRMat<Double>& coef_mat)
{
	static const Int4 ISIZE = 3;
	assert( coef_mat.nrows() == ISIZE && coef_mat.ncols() == ISIZE );
	for(Int4 k=0; k<ISIZE; k++)
	{
		const VecDat3<Double>& coef_tray = rhs[k].putCoef();
		for(Int4 k2=0; k2<ISIZE; k2++) coef_mat[k][k2] = coef_tray[k2];
	}
}

#endif /*XYZCoord_HH_*/


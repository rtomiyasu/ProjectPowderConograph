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
#include "MillerIndex.hh"

MillerIndex3::MillerIndex3()
{
	for(Int4 k=0; k<ISIZE+1; k++) m_hkl[k] = 0;
}

MillerIndex3::MillerIndex3(const MillerIndex3& rhs)
{
	for(Int4 k=0; k<ISIZE+1; k++) m_hkl[k] = rhs.m_hkl[k];
}

MillerIndex3::MillerIndex3(const VecDat3<Int4>& rhs)
{
	for(Int4 k=0; k<ISIZE; k++) m_hkl[k] = rhs[k];
	m_hkl[3]=-m_hkl[0]-m_hkl[1];	// i
}

MillerIndex3::MillerIndex3(const Int4& h, const Int4& k, const Int4& l)
{
	m_hkl[0]=h;	// h
	m_hkl[1]=k;	// k
	m_hkl[2]=l;	// l
	m_hkl[3]=-h-k;	// i
}

MillerIndex3::~MillerIndex3()
{
}

MillerIndex3& MillerIndex3::operator=(const MillerIndex3& rhs)
{
	if(this != &rhs)
		for(Int4 k=0; k<ISIZE+1; k++) m_hkl[k] = rhs.m_hkl[k];
	return *this;
}

inline Double norm(const MillerIndex3& hkl, const SymMat<Double>& S)
{
	return hkl[0]*hkl[0]*S(0,0)+hkl[1]*hkl[1]*S(1,1)+hkl[2]*hkl[2]*S(2,2)
	+ 2.0*( hkl[0]*hkl[1]*S(0,1)+hkl[0]*hkl[2]*S(0,2)+hkl[1]*hkl[2]*S(1,2) );
}

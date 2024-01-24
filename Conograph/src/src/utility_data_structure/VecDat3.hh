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
#ifndef VECDAT3_H_
#define VECDAT3_H_

#include<complex>

using namespace std;

// Class of a vector of size 3 with some operators.
template <class T>
class VecDat3
{
protected:
	enum{ ISIZE = 3 };
	T m_vec[ISIZE];
	
public:
	VecDat3();
	VecDat3(const T&);
	VecDat3(const T&, const T&, const T&);
	VecDat3(const VecDat3<T>&); // copy constructor
	virtual ~VecDat3();
	VecDat3<T>& operator=(const T&);
	VecDat3<T>& operator=(const VecDat3<T>&);
	inline VecDat3<T>& operator+=(const VecDat3<T>&);
	inline VecDat3<T>& operator-=(const VecDat3<T>&);
	inline VecDat3<T>& operator*=(const T&);
	inline VecDat3<T> operator*(const T&) const;
	inline T& operator[](const int&);
	inline const T& operator[](const int&) const;
};

template <class T>
VecDat3<T>::VecDat3()
{
//	for(int k=0; k<ISIZE; k++) m_vec[k]=0;
}

template <class T>
VecDat3<T>::VecDat3(const T& a)
{
	for(int k=0; k<ISIZE; k++) m_vec[k]=a;
}

template <class T>
VecDat3<T>::VecDat3(const T& a, const T& b, const T& c)
{
	m_vec[0]=a;
	m_vec[1]=b;
	m_vec[2]=c;
}

template <class T>
VecDat3<T>::VecDat3(const VecDat3<T> &rhs)
{
	for(int k=0; k<ISIZE; k++) m_vec[k]=rhs[k];
}

template <class T>
VecDat3<T>::~VecDat3()
{
}

template <class T>
VecDat3<T>& VecDat3<T>::operator=(const T& a)
{
	for(int k=0; k<ISIZE; k++) m_vec[k]=a;
	return *this;
}

template <class T>
VecDat3<T>& VecDat3<T>::operator=(const VecDat3<T>& rhs)
{
	if(this != &rhs)
		for(int k=0; k<ISIZE; k++) m_vec[k]=rhs[k];
	return *this;
}

template <class T>
inline VecDat3<T>& VecDat3<T>::operator+=(const VecDat3<T>& rhs)
{
	for(int k=0; k<ISIZE; k++) m_vec[k]+=rhs[k];
	return *this;
}

template <class T>
inline VecDat3<T>& VecDat3<T>::operator-=(const VecDat3<T>& rhs)
{
	for(int k=0; k<ISIZE; k++) m_vec[k]-=rhs[k];
	return *this;
}

template <class T>
inline VecDat3<T>& VecDat3<T>::operator*=(const T& a)
{
	for(int k=0; k<ISIZE; k++) m_vec[k]*=a;
	return *this;
}

template <class T>
inline VecDat3<T> VecDat3<T>::operator*(const T& a) const
{
	VecDat3<T> ans(*this);
	for(int k=0; k<ISIZE; k++) ans[k]*=a;
	return ans;
}

template <class T>
inline T& VecDat3<T>::operator[](const int& k)
{
	return m_vec[k];
}

template <class T>
inline const T& VecDat3<T>::operator[](const int& k) const
{
	return m_vec[k];
}

template <class T>
inline VecDat3<T> operator+(const VecDat3<T>& lhs, const VecDat3<T>& rhs)
{
	VecDat3<T> ans(lhs);
	ans+=rhs;
	return ans;
}

template <class T>
inline VecDat3<T> operator-(const VecDat3<T>& lhs, const VecDat3<T>& rhs)
{
	VecDat3<T> ans(lhs);
	ans-=rhs;
	return ans;
}

template <class T>
inline bool operator<(const VecDat3<T>& lhs, const VecDat3<T>& rhs)
{
	if( lhs[0] < rhs[0] ) return true;
	if( rhs[0] < lhs[0] ) return false;

	if( lhs[1] < rhs[1] ) return true;
	if( rhs[1] < lhs[1] ) return false;

	if( lhs[2] < rhs[2] ) return true;
	return false;
}


#endif /*VECDAT3_H_*/

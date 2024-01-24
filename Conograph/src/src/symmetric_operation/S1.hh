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
#ifndef S1_HH_
#define S1_HH_

#include <complex>
#include <iosfwd>
#include "../RietveldAnalysisTypes.hh"
#include "../utility_func/zmath.hh"

// The class of the complex number with absolute value 1.
class S1
{
	friend inline S1 operator-(const S1& rhs);
	friend inline bool operator==(const S1& lhs, const S1& rhs);
	friend inline S1 operator*(const S1& a, const Int4& k);
//	friend inline S1 operator*(const S1& a, const Double& k);
	friend inline CMPX_DP ipow_abs1(const S1& a, const Int4& k);

private:
	CMPX_DP mem;

public:
	S1();
	S1(const S1&);
	S1(const Double&);
	S1(const Int4&, const Int4&);
	~S1();

	S1& operator=(const S1&);
	Double toDouble() const;
	string toString() const;

	inline S1& operator+=(const S1&);
	inline S1& operator-=(const S1&);
};

inline S1& S1::operator+=(const S1& rhs)
{
	mem *= rhs.mem;
	return *this;
}

inline S1& S1::operator-=(const S1& rhs)
{
	mem *= conj(rhs.mem);
	return *this;
}

inline S1 operator+(const S1& lhs, const S1& rhs)
{
	S1 ans(lhs);
	ans += rhs;
	return ans;
}

inline S1 operator-(const S1& lhs, const S1& rhs)
{
	S1 ans(lhs);
	ans -= rhs;
	return ans;
}

inline S1 operator-(const S1& rhs)
{
	S1 ans;
	ans.mem = conj(rhs.mem);
	return ans;
}

inline S1 operator*(const S1& lhs, const Int4& rhs)
{
	S1 ans;
	if(rhs>=0) ans.mem = uipow(lhs.mem, (unsigned char)rhs);
	else ans.mem = uipow(conj(lhs.mem), (unsigned char)(-rhs));
	return ans;
}

//inline S1 operator*(const S1& lhs, const Double& rhs)
//{
//	return S1( lhs.toDouble() * rhs );
//}

//inline bool operator<(const S1& lhs, const S1& rhs)
//{
//	static const Double thred = 1.0e-2;
//	if( lhs.toDouble() + thred < rhs.toDouble() ) return true;
//	return false;
//}

inline bool operator==(const S1& lhs, const S1& rhs)
{
	static const Double thred = 1.0e-2;
	if( abs(lhs.mem - rhs.mem) > thred) return false;
	return true;
}

//inline bool operator!=(const S1& lhs, const S1& rhs)
//{
//	return !(lhs == rhs);
//}

// (The absolute value of the argument needs to equals 1.)
inline CMPX_DP ipow_abs1(const S1& a, const Int4& k)
{
	S1 ans = a * k;
	return ans.mem;
}

#endif /*S1_HH_*/


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
#ifndef StringS1_HH_
#define StringS1_HH_

#include<iosfwd>
#include"../RietveldAnalysisTypes.hh"
#include"../utility_func/zstring.hh"
#include"../utility_func/gcd.hh"
#include"../utility_data_structure/VecDat3.hh"
#include"S1.hh"

// Class of the data in the form of ax + by + cz + d.
class StringS1
{
	friend inline StringS1 operator-(const StringS1& rhs);
	friend inline bool operator<(const StringS1& lhs, const StringS1& rhs);
	friend inline bool operator==(const StringS1& lhs, const StringS1& rhs);
//	friend inline StringS1 operator*(const XYZCoord<StringS1>& lhs, const StringS1& rhs);
//	friend inline XYZCoord<StringS1> inverse(const XYZCoord<StringS1>& lhs);

private:
	VecDat3<Double> m_xyz_coef; // a, b, c
	S1 m_s1;	// d

public:
	StringS1();
	StringS1(const StringS1&);
	StringS1(const S1&);
	StringS1(const VecDat3<Double>&, const S1&);
	~StringS1();

	StringS1& operator=(const StringS1&);
	StringS1& operator=(const S1&);
	inline StringS1& operator+=(const StringS1&);
	inline StringS1 operator+(const StringS1&) const;
	inline StringS1& operator-=(const StringS1&);
	inline StringS1 operator-(const StringS1&) const;

	inline string toString() const;
	inline const S1& toS1() const;
	inline const VecDat3<Double>& putCoef() const;
	
	static const StringS1& putX();
	static const StringS1& putY();
	static const StringS1& putZ();
};

inline StringS1& StringS1::operator+=(const StringS1& rhs)
{
	m_xyz_coef += rhs.m_xyz_coef;
	m_s1 += rhs.m_s1;
	return *this;
}

inline StringS1 StringS1::operator+(const StringS1& rhs) const
{
	StringS1 ans(*this);
	ans.m_xyz_coef += rhs.m_xyz_coef;
	ans.m_s1 += rhs.m_s1;
	return ans;
}

inline StringS1& StringS1::operator-=(const StringS1& rhs)
{
	m_xyz_coef -= rhs.m_xyz_coef;
	m_s1 -= rhs.m_s1;
	return *this;
}

inline StringS1 StringS1::operator-(const StringS1& rhs) const
{
	StringS1 ans(*this);
	ans.m_xyz_coef -= rhs.m_xyz_coef;
	ans.m_s1 -= rhs.m_s1;
	return ans;
}

inline StringS1 operator*(const StringS1& lhs, const Int4& rhs)
{
	return StringS1( lhs.putCoef() * rhs, lhs.toS1() * rhs );
}

//inline StringS1 operator*(const StringS1& lhs, const Double& rhs)
//{
//	return StringS1( lhs.putCoef() * rhs, lhs.toS1() * rhs );
//}

inline string StringS1::toString() const
{
	static const string xyz[] = {"x","y","z"};
	string str[4];
	for(Int4 k=0; k<3; k++)
	{
		pair<Int4, Int4> frac;
		if( !dbl2fraction(m_xyz_coef[k], frac) ) str[k] = num2str<Double>( m_xyz_coef[k] );
		else if( frac.first == 0 ) continue;
		else if( frac.second == 1 )
		{
			if( frac.first == -1 ) str[k] = "-";
			else if( frac.first != 1 ) str[k] = num2str<Int4>( frac.first );
		}
		else
		{
			str[k] = num2str<Int4>(frac.first)+"/"+num2str<Int4>(frac.second);
		}
		str[k] += xyz[k];
	}
	str[3] = m_s1.toString();
	if( str[3] == "0" || str[3] == "1" ) str[3].clear();

	string str2 = str[0];
	for(Int4 k=1; k<4; k++){
		if( str[k].empty() ) continue;
		if( str2.empty() || str[k].at(0) == '-' || str[k].at(0) == '+' ) str2 += str[k];
		else str2 += "+" + str[k];
	}
	
	if(str2.empty()) return "0";
	else return str2;
}


inline const VecDat3<Double>& StringS1::putCoef() const
{
	return m_xyz_coef;
}

inline const S1& StringS1::toS1() const
{
	return m_s1;
}


inline StringS1 operator-(const StringS1& rhs)
{
	StringS1 ans( rhs.m_xyz_coef*(-1), -rhs.m_s1);
	return ans;
}

//inline bool operator<(const StringS1& lhs, const StringS1& rhs)
//{
//	static const Double thred = 1.0e-10;
//	for(Int4 i=0; i<3; i++)
//	{
//		if( lhs.m_xyz_coef[i] + thred < rhs.m_xyz_coef[i] ) return true;
//		if( rhs.m_xyz_coef[i] + thred < lhs.m_xyz_coef[i] ) return false;
//	}
//	if( lhs.m_s1 < rhs.m_s1 ) return true;
//	return false;
//}

inline bool operator==(const StringS1& lhs, const StringS1& rhs)
{
	static const Double thred = 1.0e-10;
	for(Int4 i=0; i<3; i++)
		if( fabs( lhs.m_xyz_coef[i] - rhs.m_xyz_coef[i] ) > thred ) return false;
	if( !(lhs.m_s1 == rhs.m_s1) ) return false;
	return true;
}

//
//inline bool operator!=(const StringS1& lhs, const StringS1& rhs)
//{
//	return !(lhs == rhs);
//}

/*
inline StringS1 operator*(const XYZCoord<StringS1>& lhs, const StringS1& rhs)
{
	VecDat3<Double> coef = 0;
	S1 c=rhs.m_s1; 
	for(Int4 k=0; k<3; k++)
	{
		for(Int4 l=0; l<3; l++) coef[l] += lhs[k].m_xyz_coef[l]*rhs.m_xyz_coef[k]; 
		c += S1( lhs[k].m_s1.toDouble() * rhs.m_xyz_coef[k] ); 
	}
		
	return StringS1(coef, c);
}


inline XYZCoord<StringS1> operator*(const XYZCoord<StringS1>& lhs, const XYZCoord<StringS1>& rhs)
{
	return XYZCoord<StringS1>(lhs*rhs[0], lhs*rhs[1], lhs*rhs[2]);
}

inline XYZCoord<StringS1> inverse(const XYZCoord<StringS1>& lhs)
{
	const Double det23 = lhs[1].m_xyz_coef[1] * lhs[2].m_xyz_coef[2] - lhs[2].m_xyz_coef[1] * lhs[1].m_xyz_coef[2];
	const Double det13_23 = lhs[1].m_xyz_coef[0] * lhs[2].m_xyz_coef[2] - lhs[2].m_xyz_coef[0] * lhs[1].m_xyz_coef[2];
	const Double det12_23 = lhs[1].m_xyz_coef[0] * lhs[2].m_xyz_coef[1] - lhs[2].m_xyz_coef[0] * lhs[1].m_xyz_coef[1];

	const Double det = 1.0 / (lhs[0].m_xyz_coef[0]*det23 - lhs[0].m_xyz_coef[1]*det13_23 + lhs[0].m_xyz_coef[2]*det12_23);

	StringS1 ans[3];
	ans[0].m_xyz_coef[0] = det23 * det; 
	ans[1].m_xyz_coef[0] = -det13_23 * det; 
	ans[2].m_xyz_coef[0] = det12_23 * det; 

	ans[0].m_xyz_coef[1] = -(lhs[0].m_xyz_coef[1] * lhs[2].m_xyz_coef[2] - lhs[2].m_xyz_coef[1] * lhs[0].m_xyz_coef[2]) * det; 
	ans[1].m_xyz_coef[1] = (lhs[0].m_xyz_coef[0] * lhs[2].m_xyz_coef[2] - lhs[2].m_xyz_coef[0] * lhs[0].m_xyz_coef[2]) * det; 
	ans[2].m_xyz_coef[1] = -(lhs[0].m_xyz_coef[0] * lhs[2].m_xyz_coef[1] - lhs[2].m_xyz_coef[0] * lhs[0].m_xyz_coef[1]) * det; 

	ans[0].m_xyz_coef[2] = (lhs[0].m_xyz_coef[1] * lhs[1].m_xyz_coef[2] - lhs[1].m_xyz_coef[1] * lhs[0].m_xyz_coef[2]) * det; 
	ans[1].m_xyz_coef[2] = -(lhs[0].m_xyz_coef[0] * lhs[1].m_xyz_coef[2] - lhs[1].m_xyz_coef[0] * lhs[0].m_xyz_coef[2]) * det; 
	ans[2].m_xyz_coef[2] = (lhs[0].m_xyz_coef[0] * lhs[1].m_xyz_coef[1] - lhs[1].m_xyz_coef[0] * lhs[0].m_xyz_coef[1]) * det; 
	
	for(Int4 i=0; i<3; i++)
	{
		ans[i].m_s1 = 0.0;
		for(Int4 k=0; k<3; k++)
		{
			ans[i].m_s1 -= S1( lhs[k].m_s1.toDouble() * ans[i].m_xyz_coef[k] ); 
		}
	}
	return XYZCoord<StringS1>(ans[0], ans[1], ans[2]);
}
*/

#endif /*StringS1_HH_*/


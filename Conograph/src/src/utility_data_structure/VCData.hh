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
#ifndef VCData_HH_
#define VCData_HH_

#ifdef DEBUG
	#include <iostream>
#endif
#include <map>
#include <assert.h>
#include "../RietveldAnalysisTypes.hh"
#include "../utility_data_structure/qdata.hh"
#include "../zerror_type/error_out.hh"
#include "map_int_int_operator.hh"

typedef Int4 type_coef;

class ZParawError;

class VCData
{
	friend inline Double calVariance(const VCData& vc);
	friend inline Double calCovariance(const VCData& vc1, const VCData& vc2);
	
private:
	enum{ m_data_num = 1 };
	enum{ m_denom = 5184 };
	static vector<QData> m_peak_qdata;
	static vector<Int4> m_peak_qindex;

	Double m_val;
	map<Int4, type_coef> m_vec_coef;

	Double calValue() const;

public:
    VCData();
    VCData(const Int4& rhs);	// If rhs = 0, set to zero.
    VCData(const Int4&, const type_coef&);
    VCData(const VCData& rhs);
    ~VCData();

	inline Double Value() const { return m_val; }; 
	inline Double Variance() const;
	inline Double Norm() const;

	VCData& operator=(const VCData&);
	inline VCData& operator+=(const VCData&);
	inline VCData& operator-=(const VCData&);
	inline VCData& operator*=(const type_coef&);
	inline VCData& operator/=(const type_coef&);

	inline const map<Int4, type_coef>& putVecCoef() const;
	inline bool operator<(const VCData&) const;
	inline bool operator<=(const VCData&) const;
	
	inline Int4 putIndex() const { if( m_vec_coef.size() == 1 ) return m_vec_coef.begin()->first; else return -1; };


	static const Int4& putDataNum();
	static Int4 putPeakNum();
	static const QData& putPeakPos(const Int4& index);

	// On input, pos_var gives the variances of pos.
	static void setQData(const vector<QData>& pos, const vector<Int4>& peak_index);
	
	static const vector<QData>& putPeakQData();
	static const vector<Int4>& putPeakQIndex();

	static type_coef putDenom() { return m_denom; }
	static VCData putVCData(const Double& rhs);
	static VCData putPeakPosArrayGuard();

	// For GUI.
	const double              &getref_m_val() const {return m_val;}
	      double              &getref_m_val()       {return m_val;}
	const map<Int4,type_coef> &getref_m_vec_coef() const {return m_vec_coef;}
	      map<Int4,type_coef> &getref_m_vec_coef()       {return m_vec_coef;}
};


inline Double calVariance(const VCData& vc)
{
	const vector<QData>& qdata = VCData::m_peak_qdata;

	Double ans=0.0;
	for(map<Int4, type_coef>::const_iterator it = vc.putVecCoef().begin(); it != vc.putVecCoef().end(); it++)
	{
		ans += qdata[ it->first ].q_var * it->second * it->second;
	}
	return ans/Double(VCData::putDenom()*VCData::putDenom());
}


inline Double calCovariance(const VCData& vc1, const VCData& vc2)
{
	const vector<QData>& qdata = VCData::m_peak_qdata;
	
	Double ans=0.0;
	map<Int4, type_coef>::const_iterator it2;
	for(map<Int4, type_coef>::const_iterator it = vc1.putVecCoef().begin(); it != vc1.putVecCoef().end(); it++)
	{
		it2 = vc2.putVecCoef().find( it->first );
		if( it2 != vc2.putVecCoef().end() ) ans += qdata[ it->first ].q_var * it->second * it2->second;
	}
	return ans/Double(VCData::putDenom()*VCData::putDenom());
}


inline Double VCData::Variance() const
{
	return calVariance(*this);
}


inline Double VCData::Norm() const
{
	Double val = Value();
	Double nor = calVariance(*this);
	return val * val / nor;
}


inline VCData& VCData::operator+=(const VCData& rhs)
{
	m_vec_coef += rhs.m_vec_coef;
	if( m_vec_coef.empty() ) m_val = 0.0;
	else m_val += rhs.m_val;
	return *this;
}

inline VCData& VCData::operator-=(const VCData& rhs)
{
	m_vec_coef -= rhs.m_vec_coef;
	if( m_vec_coef.empty() ) m_val = 0.0;
	else m_val -= rhs.m_val;
	return *this;
}

inline VCData& VCData::operator*=(const type_coef& rhs)
{
	m_val *= rhs;
	m_vec_coef *= rhs;
	return *this;
}

inline VCData& VCData::operator/=(const type_coef& rhs)
{
#ifdef DEBUG
	for(map<Int4, type_coef>::const_iterator it=m_vec_coef.begin(); it!=m_vec_coef.end(); it++)
	{
		if( it->second % rhs != 0 ) cout << it->second << " " << rhs << endl;
		assert( it->second % rhs == 0 );
	}
#endif
	
	m_val /= rhs;
	m_vec_coef /= rhs;
	return *this;
}

inline const map<Int4, type_coef>& VCData::putVecCoef() const
{
	return m_vec_coef;
}


inline VCData operator+(const VCData& lhs, const VCData& rhs)
{
	VCData ans(lhs);
	ans += rhs;
	return ans;
}

inline VCData operator-(const VCData& lhs, const VCData& rhs)
{
	VCData ans(lhs);
	ans -= rhs;
	return ans;
}

inline VCData operator*(const VCData& lhs, const type_coef& rhs)
{
	VCData ans(lhs);
	ans *= rhs;
	return ans;
}

inline VCData operator/(const VCData& lhs, const type_coef& rhs)
{
	VCData ans(lhs);
	ans /= rhs;
	return ans;
}

inline bool VCData::operator<(const VCData& rhs) const
{
//	if( this->m_vec_coef == rhs.m_vec_coef ) return false;
	return ( this->Value() < rhs.Value() );
}

inline bool VCData::operator<=(const VCData& rhs) const
{
//	if( this->m_vec_coef == rhs.m_vec_coef ) return false;
	return ( this->Value() <= rhs.Value() );
}

inline ostream& operator<<(ostream& os, const VCData& rhs)
{ 
	os << rhs.Value();
	return os;
}


inline bool vc_equiv(const VCData& lhs, const VCData& rhs, const Double& cv2)
{
	const VCData diff = lhs - rhs;
	return ( diff.Value()*diff.Value() <= cv2*min(lhs.Variance(), rhs.Variance() ) );
}

template<class Iterator>
pair<Iterator, Iterator> equal_range(const Iterator& it_begin,
		const Iterator& it_end, const VCData& rhs, const Double& cv2)
{
	const Double dval = rhs.Value();
	const Double dvar = sqrt( cv2*rhs.Variance() );
	pair<Iterator, Iterator> ans;
	ans.first = lower_bound( it_begin, it_end, dval - dvar );
	ans.second = upper_bound( ans.first, it_end, dval + dvar );
	if( ans.second < it_end ) ans.second++;
	return ans;
}


#endif /*VCData_HH_*/

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
#include"VCData.hh"
#include"../utility_func/gcd.hh"
#include"../zparam/ZParawError.hh"

using namespace std;

vector<QData> VCData::m_peak_qdata;
vector<Int4> VCData::m_peak_qindex;

Double VCData::calValue() const
{
	Double ans = 0.0;
	for(map<Int4, type_coef>::const_iterator it = m_vec_coef.begin(); it != m_vec_coef.end(); it++)
	{
		ans += m_peak_qdata[it->first].q * it->second;
	}
	return ans/Double(m_denom);
}


VCData::VCData()
{
	m_val = 0.0;
}


VCData::VCData(const Int4& rhs) // Constructor for sort functions.
{
	assert( rhs == 0 );
	m_val = rhs;
}

VCData::VCData(const Int4& peak_index, const type_coef& coef)
{
	assert( 0 <= peak_index && peak_index < (Int4)m_peak_qdata.size() );

	m_vec_coef.insert( map<Int4, type_coef>::value_type( peak_index, coef*m_denom ) );
	m_val = this->calValue();
}


VCData::VCData(const VCData& rhs)
{
	m_vec_coef = rhs.m_vec_coef;
	m_val = rhs.m_val;
}


VCData::~VCData()
{
}


VCData& VCData::operator=(const VCData& rhs)
{
	if (this != &rhs)
	{
		m_vec_coef = rhs.m_vec_coef;
		m_val = rhs.m_val;
	}
	return *this;
}


Int4 VCData::putPeakNum()
{
	return m_peak_qdata.size();
}

const QData& VCData::putPeakPos(const Int4& index)
{
	assert( 0 <= index && index < (Int4)m_peak_qdata.size() );
	return m_peak_qdata[index];
}


void VCData::setQData(const vector<QData>& arg1, const vector<Int4>& arg2)
{
	assert( arg1.size() == arg2.size() );
	m_peak_qdata = arg1;
	m_peak_qindex = arg2;
}


const vector<QData>& VCData::putPeakQData()
{
	return m_peak_qdata;
}

const vector<Int4>& VCData::putPeakQIndex()
{
	return m_peak_qindex;
}

VCData VCData::putVCData(const Double& rhs)
{
	VCData ans;
	ans.m_val = rhs;
	return ans;
}

VCData VCData::putPeakPosArrayGuard()
{
	VCData ans;
	ans.m_vec_coef.insert( map<Int4, type_coef>::value_type( m_peak_qdata.size(), m_denom ) );
	ans.m_val = 0.0;
	return ans;
}

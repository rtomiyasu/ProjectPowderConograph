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
#ifndef QDATA_HH_
#define QDATA_HH_

#include <assert.h>
#include <algorithm>
#include "../utility_func/zmath.hh"

class QData
{
public:
	Double q;
	Double q_var;

	QData(){};
	QData(const Double& arg1, const Double& arg2){ q = arg1; q_var = arg2; };

	inline bool operator<(const QData& rhs) const { return q < rhs.q; };
//	inline bool operator<(const Double& rhs) const { return q < rhs; };
};

inline Double put_minimum_counitcell_edge_length(const vector<QData>& qdata, const Int4& used_peak_num)
{
	const Int4 qnum = qdata.size();
	assert( qnum > 1 && used_peak_num > 1 );

	const Int4 CHECK_PEAK_NUM = (used_peak_num>qnum?qnum:used_peak_num);
	const Double Q_MAX = (qdata.begin()+(CHECK_PEAK_NUM-1))->q;	// We assume here Qdata is sorted into ascending order.
	const Double Q_MIN = qdata.begin()->q;

	return (sqrt(Q_MAX) - sqrt(Q_MIN)) / (CHECK_PEAK_NUM - 1);
}


inline Double put_maximum_counitcell_volume(const vector<QData>& qdata, const Int4& used_peak_num)
{
	const Int4 qnum = qdata.size();
	assert( qnum > 1 && used_peak_num > 1 );

	const Int4 CHECK_PEAK_NUM = (used_peak_num>qnum?qnum:used_peak_num);
	const Double Q_MAX = (qdata.begin()+(CHECK_PEAK_NUM-1))->q;	// We assume here Qdata is sorted into ascending order.
	const Double Q_MIN = qdata.begin()->q;
	return ( PI() * ( Q_MAX*sqrt(Q_MAX) - Q_MIN*sqrt(Q_MIN) ) ) / ( (CHECK_PEAK_NUM - 1) * 1.5);
}


inline pair<vector<QData>::const_iterator, vector<QData>::const_iterator> closest_qdata(
		const vector<QData>::const_iterator& it_begin,
		const vector<QData>::const_iterator& it_end,
		const Double& rhs, const Double& rhs_diff)
{
	const QData qrhs(rhs, 0.0);
	pair<vector<QData>::const_iterator, vector<QData>::const_iterator> it_pair_dbl2;
	it_pair_dbl2.first = lower_bound( it_begin, it_end, qrhs );
	it_pair_dbl2.second = upper_bound( it_pair_dbl2.first, it_end, qrhs );
	if( it_pair_dbl2.first >= it_pair_dbl2.second )
	{
		if( it_pair_dbl2.first > it_begin )
		{
			if( (it_pair_dbl2.first - 1)->q >= rhs - rhs_diff ) it_pair_dbl2.first--;
		}

		if( it_pair_dbl2.second < it_end )
		{
			if( it_pair_dbl2.second->q <= rhs + rhs_diff ) it_pair_dbl2.second++;
		}
	}

	return it_pair_dbl2;
}

inline pair<vector<QData>::const_iterator, vector<QData>::const_iterator> equal_range_qdata(
		const vector<QData>::const_iterator& it_begin,
		const vector<QData>::const_iterator& it_end,
		const Double& rhs, const Double& rhs_diff)
{
	pair<vector<QData>::const_iterator, vector<QData>::const_iterator> it_pair_dbl2;
	it_pair_dbl2.first = lower_bound( it_begin, it_end, QData(rhs - rhs_diff, 0.0) );
	it_pair_dbl2.second = upper_bound( it_pair_dbl2.first, it_end, QData(rhs + rhs_diff, 0.0) );
	return it_pair_dbl2;
}

#endif /* QDATA_HH_ */

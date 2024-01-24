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
#ifndef RANGE_HH_
#define RANGE_HH_

#include <algorithm>
#include <limits>

#include "../RietveldAnalysisTypes.hh"

template <class T>
class range{
public:
	T begin;
	T end;
	range(){};
	range(const T& a, const T& b){ begin=a; end=b; };
};

template <class T>
inline void put_range_union(const range<T>& lhs, const range<T>& rhs, range<T>& ans)
{
	ans.begin = min(lhs.begin,rhs.begin);
	ans.end = max(lhs.end,rhs.end);
}

template <class T>
void put_range_union(const vector< range<T> >& d_range2, range<T>& d_range)
{
	static const T MAX_DP = numeric_limits<T>::max();

	const Int4 range_num = d_range2.size();

	d_range.begin = MAX_DP;
	d_range.end = -MAX_DP;
	
	for(Int4 k=0; k<range_num; k++)
	{
		put_range_union(d_range, d_range2[k], d_range);
	}
}


#endif /*RANGE_HH_*/

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
#ifndef _INDEX_SET_HH_
#define _INDEX_SET_HH_
#include<vector>
#include<algorithm>
#include "../RietveldAnalysisTypes.hh"
#include "../zparam/etype_ID.hh"

template <class T>
class index_set{
public:
	Int4 index;
	T element;
	index_set(){};
	index_set(const Int4& i, const T& v) : index(i), element(v) {};
};

// Vec_DP_save are used as a container of vectors for the saving memory.
typedef vector< index_set<Double> > Vec_DP_save;

struct constr_DP{
	etype_ID ID;
	Vec_DP_save constr;
};

// Mat_DP_save are used as a container of flags and constraints for the saving memory.
typedef vector< constr_DP > Mat_DP_constr;


// Returns the inner product.(lhs, rhs are vectors.)
inline Double product(const Vec_DP_save& lhs, const Double* rhs)
{
	Double t = 0.0;
	for(UInt4 l=0; l<lhs.size(); l++) t += rhs[lhs[l].index] * lhs[l].element;
	return t;
}

// Change the type of a vector from Vec_DP_save to Vec_DP.
inline void chToVector(const Vec_DP_save& lhs, const Int4 pn, Vec_DP& rhs)
{
 	rhs.clear();
 	rhs.resize(pn, 0.0);
 	for(UInt4 j=0; j<lhs.size(); j++) rhs[lhs[j].index] += lhs[j].element;
}


template <class T>
bool operator<(const index_set<T>& lhs, const index_set<T>& rhs)
{
	return lhs.element < rhs.element;
}

template <class T>
bool cmp_index(const index_set<T>& lhs, const index_set<T>& rhs)
{
	return lhs.index < rhs.index;
}


// Permute the elements. 
template<class Iterator>
void permute(const Vec_INT& order, const Iterator& it_start)
{
    const int isize = order.size();

	Vec_BOOL vflag(isize+1, true);	// The last index is against array overflow.
	
	Int4 row_start = 0, i, i2;
	while(row_start < isize)
	{
		i = row_start++;
		i2 = order[i];
		vflag[i] = false;
		if(i != i2)
		{
			while( vflag[i2] )
			{
				iter_swap( it_start + i, it_start + i2 );
				vflag[i2] = false;
				i = i2;
				i2 = order[i];
			}
		}

		while( !vflag[row_start] ) row_start++;
	}
}



template <class S>
void sort_order(vector<S>& lhs, Vec_INT& order)
{
	const Int4 isize = lhs.size();
	
	vector< index_set<S> > lhstray;
	for(Int4 k=0; k<isize; k++)
	{
		lhstray.push_back( index_set<S>(k, lhs[k]) );
	}

	stable_sort(lhstray.begin(), lhstray.end());

	order.resize(isize);
	for(Int4 k=0; k<isize; k++)
	{
		order[k] = lhstray[k].index;
		lhs[k] = lhstray[k].element;
	}
}

#endif


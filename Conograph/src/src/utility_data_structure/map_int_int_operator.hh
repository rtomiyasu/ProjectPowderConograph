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
#ifndef MAP_INT_INT_OPERATOR_HH_
#define MAP_INT_INT_OPERATOR_HH_

#include<map>
#include<algorithm>
#include"../RietveldAnalysisTypes.hh"
#include"../utility_func/zstring.hh"

template<class T>
inline map<Int4, T>& operator+=(map<Int4, T>& lhs, const map<Int4, T>& rhs)
{
	typename map<Int4, T>::const_iterator it = rhs.begin();
	typename map<Int4, T>::iterator it2;
	while( it != rhs.end() )
	{
		it2 = lhs.find( it->first );
		if( it2 != lhs.end() )
		{
			it2->second += it->second;
			if( it2->second == 0 ) lhs.erase(it2);
		}
		else lhs[it->first] = it->second;
		++it;
	}	
	return lhs;
}

template<class T>
inline map<Int4, T> operator+(const map<Int4, T>& lhs, const map<Int4, T>& rhs)
{
	map<Int4, T> ans = lhs;
	ans += rhs;
	return ans;
}

template<class T>
inline map<Int4, T>& operator-=(map<Int4, T>& lhs, const map<Int4, T>& rhs)
{
	typename map<Int4, T>::const_iterator it = rhs.begin();
	typename map<Int4, T>::iterator it2;
	while( it != rhs.end() )
	{
		it2 = lhs.find( it->first );
		if( it2 != lhs.end() )
		{
			it2->second -= it->second;
			if( it2->second == 0 ) lhs.erase(it2);
		}
		else lhs[it->first] = - it->second;
		++it;
	}	
	return lhs;
}

template<class T>
inline map<Int4, T> operator-(const map<Int4, T>& lhs, const map<Int4, T>& rhs)
{
	map<Int4, T> ans = lhs;
	ans -= rhs;
	return ans;
}

template<class T>
inline map<Int4, T>& operator*=(map<Int4, T>& lhs, const T& rhs)
{
	if (rhs == 0) lhs.clear();
	else
	{
		typename map<Int4, T>::iterator it = lhs.begin();
		while( it != lhs.end() )
		{
			it->second *= rhs;
			++it;
		}
	}
	return lhs;
}

template<class T>
inline map<Int4, T> operator*(const map<Int4, T>& lhs, const T& rhs)
{
	map<Int4, T> ans = lhs;
	ans *= rhs;
	return ans;
}

template<class T>
inline map<Int4, T>& operator/=(map<Int4, T>& lhs, const T& rhs)
{
	typename map<Int4, T>::iterator it = lhs.begin();
	while( it != lhs.end() )
	{
		it->second /= rhs;
		++it;
	}
	return lhs;
}


template<class T>
inline map<Int4, T> operator/(map<Int4, T>& lhs, const T& rhs)
{
	map<Int4, T> ans = lhs;
	ans /= rhs;
	return ans;
}


template<class T>
T inner_product(const map<Int4, T>& lhs, const map<Int4, T>& rhs)
{
	T ans = 0;
	
	typename map<Int4, T>::const_iterator it;
	typename map<Int4, T>::const_iterator it2;

	it = lhs.begin();
	while( it != lhs.end() )
	{
		it2 = rhs.find( it->first );
		if( it2 != rhs.end() ) ans += it2->second * it->second;
		++it;
	}

	return ans;
}


template<class T>
map<Int4, T> product(const map<Int4, T>& lhs, const vector< map<Int4, T> >& rhs)
{
	map<Int4, T> ans;
	T num;

	const Int4 num_col = rhs.size();
	for(Int4 k=0; k<num_col; k++)
	{
		num = inner_product(lhs, rhs[k]);
		if( num != 0 )
			ans.insert( typename map<Int4, T>::value_type( k, num ) );
	}	

	return ans;
}



// void chToVector(const map<Int4, T>&, const Int4, Vec_DP&);

template<class T>
string chToString(const map<Int4, T>& lhs, const string& str)
{
 	typename map<Int4, T>::const_iterator it = lhs.begin();
 	if( it == lhs.end() ) return "0";
 	
	string ans;
	if( it->second > 0 )
	{
		if( it->second!=1 ) ans+= num2str<T>( it->second );
	}
	else
	{
		if( it->second==-1 ) ans += "-";
		else ans += num2str<T>( it->second );
	}
	ans += str+"_"+num2str<T>( it->first+1 );
	it++;
 	
 	while( it != lhs.end() )
 	{
 		if( it->second > 0 )
 		{
			ans += "+";
 			if( it->second!=1 ) ans+= num2str<T>( it->second );
 		}
 		else
 		{
 			if( it->second==-1 ) ans += "-";
 			else ans += num2str<T>( it->second );
 		}
 		ans += str+"_"+num2str<T>( it->first+1 );
 		it++;
 	}
 	return ans;
}


template<class T>
T sum(const typename map<Int4, T>::const_iterator& it_begin, const typename map<Int4, T>::const_iterator& it_end)
{
	T ans = 0;
	typename map<Int4, T>::const_iterator it = it_begin;
 	while( it != it_end )
 	{
		ans += it->second;
		it++;
 	}
 	return ans;
}


template<class T>
bool operator==(const map<Int4,T>& lhs, const map<Int4,T>& rhs)
{
	if( lhs.size() != rhs.size() ) return false;
	const map<Int4,T> diff = rhs - lhs;
	return diff.empty();
}


template<class T>
bool cmp(const map<Int4,T>& lhs, const map<Int4,T>& rhs)
{
	if( lhs.size() < rhs.size() ) return true;
	if( lhs.size() > rhs.size() ) return false;
	const map<Int4,T> diff = rhs - lhs;
	if( diff.empty() ) return false;
	assert( diff.begin()->second != 0 );
	return ( diff.begin()->second > 0 );
}

#endif /*MAP_INT_INT_OPERATOR_HH_*/

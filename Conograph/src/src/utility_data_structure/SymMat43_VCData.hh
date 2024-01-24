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
#ifndef SYMMAT43_VCDATA_HH_
#define SYMMAT43_VCDATA_HH_

#include "../utility_data_structure/SymMat.hh"
#include "../utility_data_structure/VCData.hh"
#include "../utility_data_structure/nrutil_nr.hh"

// Assume the first element is a symmetric matrix of size 3 and the second element is a 4-by-3 matrix.
typedef pair< SymMat<VCData>, NRMat<Int4> > SymMat43_VCData;

inline bool operator<(const SymMat<VCData>& lhs, const SymMat<VCData>& rhs)
{
	assert( lhs.size() == rhs.size() );

	const Int4 ISIZE = lhs.size();

	for(Int4 i=0; i<ISIZE; i++)
		for(Int4 j=i; j<ISIZE; j++)
		{
			const map<Int4,type_coef> lhs_sym = lhs(i,j).putVecCoef();
			const map<Int4,type_coef> rhs_sym = rhs(i,j).putVecCoef();

			if( lhs_sym.size() < rhs_sym.size() ) return true;
			if( lhs_sym.size() > rhs_sym.size() ) return false;

			const map<Int4,type_coef> diff = rhs_sym - lhs_sym;
			if( diff.empty() ) continue;
			assert( diff.begin()->second != 0 );
			if( diff.begin()->second > 0 ) return true;
			else return false;
		}
	return false;
};


inline bool operator==(const SymMat<VCData>& lhs, const SymMat<VCData>& rhs)
{
	assert( lhs.size() >= 3 );
	assert( rhs.size() >= 3 );

	map<Int4,type_coef> diff;

	for(Int4 i=0; i<3; i++)
	{
		for(Int4 j=i; j<3; j++)
		{
			if( lhs(i,j).putVecCoef() == rhs(i,j).putVecCoef() ) continue;
			return false;
		}
	}
	return true;
};


#endif /*SYMMAT43_HH_*/

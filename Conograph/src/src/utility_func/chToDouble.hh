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
#ifndef CHTODOUBLE_HH_
#define CHTODOUBLE_HH_

#include"../utility_data_structure/SymMatWCovar.hh"
#include"../utility_data_structure/VCData.hh"

inline SymMat<Double> chToDouble(const SymMat<Double>& arg)
{
	return arg;
};

inline SymMat<Double> chToDouble(const SymMat<VCData>& vc_S)
{
	const Int4 ISIZE = vc_S.size();
	SymMat<Double> ans(ISIZE);

	for(Int4 i=0; i<ISIZE; i++)
		for(Int4 j=i; j<ISIZE; j++)
			ans(i,j) = vc_S(i,j).Value();
	return ans;
};


inline SymMatWCovar chToDoubleWCovar(const SymMat<VCData>& vc_S)
{
	const Int4 ISIZE = vc_S.size();
	
	SymMatWCovar ans(ISIZE);
	for(Int4 i=0; i<ISIZE; i++)
		for(Int4 j=i; j<ISIZE; j++)
			ans(i,j) = vc_S(i,j).Value();

	for(Int4 i=0; i<ISIZE; i++)
		for(Int4 i2=i; i2<ISIZE; i2++)
		{
			for(Int4 j=0; j<ISIZE; j++)
				for(Int4 j2=j; j2<ISIZE; j2++)
				{
					ans.Covar(i,i2,j,j2) = calCovariance(vc_S(i,i2), vc_S(j,j2));
				}
		}
	return ans;
};


#endif /*CHTODOUBLE_HH_*/

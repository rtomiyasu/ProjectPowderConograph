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
#ifndef _GCD_HH_
#define _GCD_HH_

#include"../RietveldAnalysisTypes.hh"

// Returns (m,n).
// On output, mt - ns = (m,n).
Int4 gcd(const Int4& m, const Int4& n, Int4& s, Int4& t);

Int4 gcd(const Int4&, const Int4&);

Int4 gcd(const vector<Int4>&, vector<Int4>&);

inline Int4 iround_half_up(const double& value)
{
	if (value >= 0.0)
   		return Int4(floor(+value + 0.5)+0.01);
	else
		return - Int4(floor(-value + 0.5)+0.01);
}

inline bool dbl2fraction(const Double& dbl, pair<Int4, Int4>& frac)
{
	const Int4 k = iround_half_up(dbl*48);
	if( fabs( k - dbl*48 ) >= 1.0e-8 ) return false;
	
	if( k == 0 )
	{
		frac.first = 0;
		frac.second = 1;
	}
	else
	{
		const Int4 common_divisor = gcd(k, 48);
		frac.first = k/common_divisor;
		frac.second = 48/common_divisor;
	}
	return true;
}

#endif /*GCD_HH_*/

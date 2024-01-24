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
#ifndef _ZMATH_HH_
#define _ZMATH_HH_

#include <assert.h>
#include "../RietveldAnalysisTypes.hh"
#include "../zparam/ZParawError.hh"

// Gets the solution of the polynomial by Newton-Raphson method. 
bool cal_polynomial_solution(const Vec_DP& poly_coef,	// Coefficients of the polynomial f(x).
								const Double& goal_Y,	// find x s.t. goal_Y = f(x).
								const Double& init_X,	// The initial value of x
								Double& X	// Answer.
							);

// Returns the value of the first argument raised to the power of the second argument.
template <class T>
inline T uipow(const T& a, const UInt4& k)
{
    T p=1;
    for (unsigned char bit = 0x80; bit > 0; bit >>= 1) {
        p *= p;
        if (k & bit) p *= a;
    }
    return p;
}

// Returns the value of the first argument raised to the power of the second argument.
template <class T>
inline T ipow(const T& a, const Int4& k)
{
	if(k>=0) return uipow(a, (unsigned char)k);
	else return uipow(1/a, (unsigned char)(-k));
}

inline const Double& PI()
{
	static const Double PI = 4.0*atan(1.0);
	return PI;
}

inline const Double& PI2()
{
	static const Double PI2 = 2.0 * PI();  // =2 pi.
	return PI2;
}

inline const CMPX_DP& PI2i()
{
	static const CMPX_DP PI2i(0, 2.0 * PI());  // =2 pi i
	return PI2i;
}

// The units of peak_shift_param are degree.
inline Double cal_theta2_deg(const vector<ZParawError>& peak_shift_param_rad, const Double& wlength, const Double& inv_d)
{
	static const Double DegRad = 180.0/PI();
	assert( peak_shift_param_rad.size() == 1 );
	const Double w_d2 = 0.5*wlength*inv_d;
	return ( 2.0*asin(w_d2) - peak_shift_param_rad[0].value )*DegRad;
}

// Return the polynomial value.
inline Double put_polynomial_value(const Vec_DP& poly_coef, const Double& x)
{
	Int4 k = poly_coef.size()-1;
	Double value = poly_coef[k];
	for(; k>0;) value = value*x + poly_coef[--k];
	return value;
}

inline Int8 iceil(const Double& num)
{
	return Int8( ceil(num) + 0.01 );
}

inline Int8 ifloor(const Double& num)
{
	return Int8( floor(num) + 0.01 );
}

#endif /*MATH_HH_*/

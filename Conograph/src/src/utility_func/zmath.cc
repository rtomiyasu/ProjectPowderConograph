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
#include "zmath.hh"


// Get the solution of the polynomial by Newton-Raphson method. 
bool cal_polynomial_solution(const Vec_DP& poly_coef, 
		const Double& goal_Y, const Double& init_X, Double& X)
{
	static const Int4 MaxIt = 1000;
	static const Double thred = 1.0e-12;

	const Int4 poly_order = poly_coef.size()-1;
	Vec_DP dcoef( poly_order );
	for(Int4 k=0; k<poly_order; k++) dcoef[k] = (k+1)*poly_coef[k+1];
	
	Double dX;
	X = init_X;
	Int4 k;
	for(k=0; k<MaxIt; k++){
		dX = ( put_polynomial_value(poly_coef, X) - goal_Y )/put_polynomial_value(dcoef, X);
		X -= dX;
		if( fabs(dX) <= thred*X ) break;
	}

	if( k >= MaxIt ) return false;
	return true;
}

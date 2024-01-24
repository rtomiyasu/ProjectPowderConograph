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
#include "gcd.hh"
#include "../zerror_type/error_out.hh"

// Returns (m,n).
// On output, mt - ns = (m,n).
Int4 gcd(const Int4& m, const Int4& n, Int4& s, Int4& t)
{
	Int4 ans;
	if( m < 0 )
	{
		ans = gcd( -m, n, s, t);
		t *= -1;
		return ans;
	}
	
	if( n < 0 )
	{
		ans = gcd( m, -n, s, t);
		s *= -1;
		return ans;
	}

	if( m == 0 )
	{
		if(n == 0) s = 0;
		else s = -1;
		t = 0;
		return n;
	}
	if( n == 0 )
	{
		s = 0; t = 1;
		return m;
	}

	if( m <= n )
	{
		Int4 d=n/m, r=n-m*d;
		ans = gcd(r, m, s, t);
		s += d * t;
		swap(s,t);
		s *= -1; t *= -1;
	}
	else
	{
		ans = gcd(n, m, t, s);
		s *= -1; t *= -1;
	}
//if(m*t-n*s!=ans)
//{
//	cerr << m << " ";
//	cerr << n << " ";
//	cerr << s << " ";
//	cerr << t << " ";
//	cerr << ans << endl;
//	cerr << "!!" << endl;
//}
	return ans;
} 


Int4 gcd(const Int4& lhs, const Int4& rhs)
{
	if(lhs==1 || rhs==1) return 1;
	if(lhs==0) return rhs;
	if(rhs==0) return lhs;
	if(lhs < 0)
	{
		return gcd(-lhs, rhs);
	}
	if(rhs < 0)
	{
		return gcd(lhs, -rhs);
	}
	
	if( lhs < rhs ) return gcd( lhs, rhs % lhs );
	else return gcd( lhs % rhs, rhs );
} 


Int4 gcd(const vector<Int4>& num_tray, vector<Int4>& coef)
{
	coef.clear();
	if( num_tray.empty() ) return 0;
	
	Int4 isize = num_tray.size();
	coef.resize(isize);
	if( isize == 1 )
	{
		if( num_tray[0] >= 0 )
		{
			coef[0] = 1;
			return num_tray[0];
		}
		else
		{
			coef[0] = -1;
			return -num_tray[0];
		}
	}

	Int4 ans = num_tray[0], s, t, index=1;
	coef[0] = 1;
	while( index < isize )
	{
		ans = gcd(ans, num_tray[index], s, t);	// ans = ans * t - num_tray[index] * s. 
		for(Int4 k=0; k<index; k++) coef[k] *= t;
		coef[index] = -s;
		index++;
	}

	return ans;
}

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
#include "gather_q_of_3D_lattice.hh"
#include "gather_q_of_Ndim_lattice.hh"

void associateQobsWithQcal(
		const vector<QData>::const_iterator it_begin, const vector<QData>::const_iterator it_end,
		const vector<HKL_Q>& hkl_q_tray,
		vector< vector<HKL_Q>::const_iterator >& closest_hkl_q_tray)
{
	closest_hkl_q_tray.clear();
	
	for(vector<QData>::const_iterator it=it_begin; it<it_end; it++)
	{
		closest_hkl_q_tray.push_back( closest_data(hkl_q_tray.begin(), hkl_q_tray.end(), it->q) );
	}
}


inline vector<QData>::const_iterator closest_data(
		const vector<QData>::const_iterator& it_begin,
		const vector<QData>::const_iterator& it_end,
		const Double& rhs)
{
	const vector<QData>::const_iterator it = lower_bound( it_begin, it_end, QData(rhs, 0.0) );
	if( it == it_begin ) return it;
	else if( it == it_end ) return it - 1;
	else
	{
		const Double diff1 = rhs - (it-1)->q;
		const Double diff2 = (it)->q - rhs;
		if( diff1 > diff2 ) return it;
		else return it - 1;
	}
}


void associateQcalWithQobs(
		const vector<HKL_Q>::const_iterator it_begin, const vector<HKL_Q>::const_iterator it_end,
		const vector<QData>& qobs_tray,
		vector< vector<QData>::const_iterator >& closest_qobs_tray)
{
	closest_qobs_tray.clear();

	for(vector<HKL_Q>::const_iterator it=it_begin; it<it_end; it++)
	{
		closest_qobs_tray.push_back( closest_data(qobs_tray.begin(), qobs_tray.end(), it->Q()) );
	}
}

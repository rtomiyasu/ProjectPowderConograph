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
#ifndef _gather_q_of_Ndim_lattice_HH_
#define _gather_q_of_Ndim_lattice_HH_
// set_additonal_Q.hh

#include "../lattice_symmetry/HKL_Q.hh"
#include "../utility_data_structure/SymMat.hh"

using namespace std;

void gatherQcal(const SymMat<Double>& S_super,
					const Double& maxQ,
					vector<HKL_Q>& hkl_q_tray
				);

void gatherQcal(const SymMat<Double>& S_super,
					const Double& maxQ,
					const NRMat<Int4>& transform_hkl,
					vector<HKL_Q>& hkl_q_tray
				);

inline void gatherQcal(const SymMat<Double>& S_super,
					const Double& maxQ,
					vector<Double>& qcal_tray
				)
{
	qcal_tray.clear();

	vector<HKL_Q> hkl_q_tray;
	gatherQcal(S_super, maxQ, hkl_q_tray);
	if( hkl_q_tray.empty() ) return;
	qcal_tray.push_back(hkl_q_tray.begin()->Q());
	for(vector<HKL_Q>::const_iterator it=hkl_q_tray.begin()+1; it!=hkl_q_tray.end(); it++)
	{
		if( *(qcal_tray.rbegin()) < it->Q() )
		{
			qcal_tray.push_back(it->Q());
		}
	}
}

// Assumes that the entries between it_begin and it_end are sorted.
inline vector<Double>::const_iterator closest_data(
		const vector<Double>::const_iterator& it_begin,
		const vector<Double>::const_iterator& it_end,
		const Double& rhs)
{
	const vector<Double>::const_iterator it = lower_bound( it_begin, it_end, rhs );
	if( it == it_begin ) return it;
	else if( it == it_end ) return it - 1;
	else
	{
		const Double diff1 = rhs - *(it-1);
		const Double diff2 = *it - rhs;
		if( diff1 > diff2 ) return it;
		else return it - 1;
	}
}


inline vector<HKL_Q>::const_iterator closest_data(
		const vector<HKL_Q>::const_iterator& it_begin,
		const vector<HKL_Q>::const_iterator& it_end,
		const Double& rhs)
{
	const vector<HKL_Q>::const_iterator it = lower_bound( it_begin, it_end, HKL_Q(NRVec<Int4>(), rhs) );
	if( it == it_begin ) return it;
	else if( it == it_end ) return it - 1;
	else
	{
		const Double diff1 = rhs - (it-1)->Q();
		const Double diff2 = it->Q() - rhs;
		if( diff1 > diff2 ) return it;
		else return it - 1;
	}
}


bool associateQcalWithQobs(const vector<HKL_Q>::const_iterator& it_begin,
	    					const vector<HKL_Q>::const_iterator& it_end,
							const Int4& scale_of_qcal,
							const vector<Double>& qobs_tray,
							const Double& resol);

void associateQobsWithQcal(const vector<Double>::const_iterator it_begin,
							const vector<Double>::const_iterator it_end,
							const vector<HKL_Q>& qcal_tray,
							vector< vector<HKL_Q>::const_iterator >& closest_qcal_tray
						);

vector<Double>::const_iterator associateQobsWithQcal(
							const vector<Double>::const_iterator& it_begin,
						    const vector<Double>::const_iterator& it_end,
						    const vector<HKL_Q>& qcal_tray,
						    const Double& resol,
						    vector< vector<HKL_Q>::const_iterator >& closest_qcal_tray);


#endif

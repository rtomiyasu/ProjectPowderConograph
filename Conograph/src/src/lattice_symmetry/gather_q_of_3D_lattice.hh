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
#ifndef _gather_q_of_3D_lattice_HH_
#define _gather_q_of_3D_lattice_HH_
// set_additonal_Q.hh

#include "HKL_Q.hh"
#include "../utility_data_structure/qdata.hh"
#include "../utility_data_structure/SymMat.hh"

using namespace std;

void associateQobsWithQcal(const vector<QData>::const_iterator it_begin,
							const vector<QData>::const_iterator it_end,
							const vector<HKL_Q>& hkl_tray,
							vector< vector<HKL_Q>::const_iterator >& closest_rep_tray
						);

void associateQcalWithQobs(
		const vector<HKL_Q>::const_iterator it_begin, const vector<HKL_Q>::const_iterator it_end,
		const vector<QData>& q_tray,
		vector< vector<QData>::const_iterator >& closest_q_tray);

#endif

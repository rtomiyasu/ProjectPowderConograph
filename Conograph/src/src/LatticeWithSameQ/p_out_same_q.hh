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
#ifndef _P_OUT_SAME_Q_HH_
#define _P_OUT_SAME_Q_HH_

#include <fstream>
#include "../RietveldAnalysisTypes.hh"

class ControlParam;
class LatticeMetricTensor;

using namespace std;

void printSolutions(const vector<LatticeFigureOfMerit>& lattice_result,
					const ControlParam& cdata,
					ostream * const os, const size_t number_output = 15);

// Output q-values in fname.
// If selected_lattice is not NULL, the information of selected_lattice is also output.
inline void printSolutions(const vector<LatticeFigureOfMerit>& lattice_result,
					const ControlParam& cdata,
					const string& fname, const size_t number_output = 15)
{
	ofstream ofs(fname.c_str());
    printSolutions(lattice_result, cdata, &ofs);
}

#endif

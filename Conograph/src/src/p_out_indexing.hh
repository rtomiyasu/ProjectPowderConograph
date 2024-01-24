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
#ifndef _P_OUT_INDEXING_HH_
#define _P_OUT_INDEXING_HH_

#include <fstream>
#include "RietveldAnalysisTypes.hh"
#include "utility_data_structure/index_set.hh"
#include "utility_data_structure/SymMat.hh"
#include "bravais_type/enumBravaisType.hh"
#include "lattice_symmetry/enumSortCriterion.hh"

class Bud;
class ControlParam;
class VCLatticeFigureOfMeritToCheckSymmetry;
class LatticeFigureOfMeritToDisplay;
class LatticeFigureOfMerit;
class Node3;
class OutputInfo;
class PeakPosData;
class QData;
class VCData;

using namespace std;

void printQdata(const vector<QData>& Qdata, ostream* os);

// Output q-values in fname.
inline void printQdata(const vector<QData>& Qdata, const string& fname)
{
    ofstream ofs(fname.c_str());
    printQdata(Qdata, &ofs);
    ofs.close();
}


void print_bud_data(const vector<Bud>& bud_basket, 
				const Vec_BOOL& root_flag, 
				const vector< vector<Int4> >& leftbr_tray,
				const vector< vector<Int4> >& rightbr_tray,
				const vector< vector<Int4> >& root_leftbr,
				const vector< vector<Int4> >& root_rightbr,
				const string& fname);

void print_bud_data(const vector<Bud>& bud_basket, const string& fname);

void print_node_data(const vector<Node3>& node_basket, const string& fname);

// If selected_lattice is not NULL, the information of selected_lattice is also output.
void printHKLdata(const vector<VCLatticeFigureOfMeritToCheckSymmetry> lattice_result[],
					const OutputInfo outinfo[],
					const eSortCriterion& sort_criterion,
					const ControlParam& cdata,
					const PeakPosData& pdata,
					const string& fname);

// If selected_lattice is not NULL, the information of selected_lattice is also output.
void printHKLdata(const vector<LatticeFigureOfMeritToDisplay> lattice_result,
					const ControlParam& cdata,
					const PeakPosData& pdata,
					const string& fname);

// Print an IGOR file.
void printPeakPosition(
		const ControlParam& cdata,
		const PeakPosData& pdata,
		const LatticeFigureOfMeritToDisplay& latfit,
		const vector<LatticeFigureOfMerit>& lattices_same_q,
		const string& fname);

// Print an IGOR file.
void printPeakPosition(
		const ControlParam& cdata,
		const PeakPosData& pdata,
		const vector<LatticeFigureOfMeritToDisplay>& latfit,
		const string& fname);

#endif

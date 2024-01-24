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
#ifndef _IndexingLattice_hh_
#define _IndexingLattice_hh_
// IndexingLattice.hh
#include "utility_data_structure/Bud.hh"
#include "utility_data_structure/SymMat43_VCData.hh"
#include "RietveldAnalysisTypes.hh"

class ControlParam;
class Bud2;
class Node3;
class PeakPosData;
class TreeLattice;
class VCData;
class VCLatticeFigureOfMeritToCheckSymmetry;

using namespace std;


class IndexingLattice
{
private:
   	Int4 m_max_peak_num;
   	Int4 m_max_edge_num;
   	Int4 m_max_node_num;

	Double m_cv2;
    Double m_min_unitcell_volume;
    Double m_max_unitcell_volume;
   	Double m_min_primitive_cell_edge;

    enum{ m_num_ref_figure_of_merit = 20 };
    ePeakShiftFunctionType m_etype_peak_shift;
    Double m_WlengthX;
    vector<ZParawError> m_peak_shift_param_rad;

	eConographAnalysisMode m_search_level;

	// Result
    vector<Bud> m_selected_bud_tray;
    vector< SymMat43_VCData > m_S_super;

public:
    IndexingLattice();
    ~IndexingLattice();

    inline const vector<Bud>& putSelectedBuds() const { return m_selected_bud_tray; };
    inline const vector<SymMat43_VCData>& putThreeDimTopographNodes() const { return m_S_super; };
    
    void setParam(const ControlParam&);

    // On input, remove the file extension from fname0.
	void determineTwoDimLattices(const PeakPosData& pdata,
								const string& fname0);

    // On input, remove the file extension from fname0.
	void determineThreeDimLattices(const string& fname0);
};

#endif

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
#ifndef _INDEXING_FUNC_DIM3_HH_
#define _INDEXING_FUNC_DIM3_HH_
// indexing_func_dim3.hh

#include <set>
#include "RietveldAnalysisTypes.hh"
#include "utility_data_structure/VecDat3.hh"
#include "utility_data_structure/SymMat43_VCData.hh"
#include "model_function/profile_function/global_function/enumPeakShiftFunctionType.hh"

using namespace std;

class Bud2;
class ZParawError;

void combinate_bud(const vector<Bud2>& bud_basket,
					const Double& cv2,
					const Double& min_unitcell_volume,
					const Double& max_unitcell_volume,
					const Double& min_length,
				    const Int4& num_ref_figure_of_merit,
				    const Double& MIN_WuM,
				    const ePeakShiftFunctionType& etype_peak_shift,
				    const Double& WlengthX,
				    const vector<ZParawError>& peak_shift_param_rad,
					const Int4& max_node_num,
					const eConographAnalysisMode& SearchLevel,
					set< pair<Double, SymMat43_VCData> >& node_basket_vol,
					set< pair<Double, SymMat43_VCData> >& node_basket_FOM
				);


#endif

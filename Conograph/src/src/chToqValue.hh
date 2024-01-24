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
#ifndef _SETVCData_HH_
#define _SETVCData_HH_

#include "RietveldAnalysisTypes.hh"
#include "ControlParam.hh"
#include "PeakPosData.hh"
#include "utility_data_structure/qdata.hh"

ZErrorMessage chXDataToQData_AD(const Double& pos, const Double& pos_sterr,
		const Double& peak_shift_param_rad, const Double& wlength, QData& qdata);

// On output, Qdata is sorted into ascending order.
ZErrorMessage chToqValue(const ControlParam& cData, const PeakPosData& pData,
		vector<QData>& pos_qdata, vector<Int4>& pos_qindex);

void fitZeroPointShift(const PeakPosData& pdata,
						const Int4& MAX_NUM_ITER,
						vector< pair<Int4, Int4> >& pair_of_q_index,
						vector< Double >& zero_point_shift_deg);

#endif /*SETVCData_HH_*/

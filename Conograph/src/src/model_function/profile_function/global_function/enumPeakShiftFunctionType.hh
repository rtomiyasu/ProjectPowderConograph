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
#ifndef ENUMPeakShiftFunctionType_HH_
#define ENUMPeakShiftFunctionType_HH_

#include "../../../zparam/etype_ID.hh"

typedef enum _PeakShiftFunctionType
{
	kPeakShiftFunction_Type0 = 0,
	kPeakShiftFunction_Type1 = 1,
} ePeakShiftFunctionType;


inline Int4 put_param_num(const ePeakShiftFunctionType& arg)
{
	if( arg == kPeakShiftFunction_Type0 ) return 0;
	else return 1;
}

inline void put_param_label(const ePeakShiftFunctionType& arg,
		vector<string>& param_label)
{
	param_label.resize(put_param_num(arg));
	if( arg == kPeakShiftFunction_Type0 ) return;
	else
	{
		param_label[0] = "delta_2theta";
	}
}

#endif /*ENUMPeakShiftFunctionType_HH_*/

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
#ifndef enumCentringType_HH_
#define enumCentringType_HH_

#include"../RietveldAnalysisTypes.hh"

enum eCentringType{ Prim=0, BaseX=1, BaseY=2, BaseZ=3, Face=4, Inner=5, Rhom_hex=6 };

enum eBASEaxis{ BaseA_Axis=0, BaseB_Axis=1, BaseC_Axis=2 };

inline string put_centring_name(const eCentringType& i)
{
	static const size_t NUM_C = 7;
	static const string name[NUM_C] = { "Primitive", "Base(A)", "Base(B)",
											"Base(C)", "Face",
											"Inner", "Rhombohedral(hexagonal-axis)" };

	return name[(size_t)i];
}


#endif /*enumCentringType_H_*/
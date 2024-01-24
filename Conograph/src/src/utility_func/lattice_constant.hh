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
#ifndef _LATTICE_CONSTANT_HH_
#define _LATTICE_CONSTANT_HH_

#include"../RietveldAnalysisTypes.hh"
#include"../utility_data_structure/SymMat.hh"
#include"../utility_data_structure/VecDat3.hh"
#include"../utility_data_structure/SymMatWCovar.hh"


inline const Int4& PNLatConst()
{
	static const Int4 PNLatConst = 6;	// The number of parameters for the isotropic Debye-Waller factor.
	return PNLatConst;
}

inline const Int4& PNXYZ()
{
	static const Int4 PNXYZ = 3;
	return PNXYZ;
}

// calculate S(i.e. A*,B*,C*,D*,E*,F*) from the lattice_constants a, b, c, alpha, beta, gamma(deg).
void calCoParameter(const VecDat3<Double>&, const VecDat3<Double>&, SymMat<Double>&);

// calculate the lattice_constants a, b, c, alpha, beta, gamma(deg) from S(i.e. A*,B*,C*,D*,E*,F*).
void calLatticeConstant(const SymMat<Double>&, VecDat3<Double>&, VecDat3<Double>&);

// calculate the lattice_constants a, b, c, alpha, beta, gamma(deg) from S(i.e. A*,B*,C*,D*,E*,F*).
void calLatticeConstant(const SymMatWCovar&, VecDat3<Double>&, VecDat3<Double>&, SymMat<Double>& LatConst_covar);


#endif /*LATTICE_CONSTANT_HH_*/

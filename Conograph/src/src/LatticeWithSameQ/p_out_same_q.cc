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
#include "../lattice_symmetry/LatticeFigureOfMeritToCheckSymmetry.hh"
#include "../ControlParam.hh"
#include "p_out_same_q.hh"




void printSolutions(const vector<LatticeFigureOfMerit>& lattice_result,
		const ControlParam& cdata,
		ostream* const os, const size_t number_output)
{
   	if( number_output < 0 ) return;

   	os->precision(6);
   	
   	VecDat3<Double> length_axis, angle_axis;
	Int4 count = 0;
	for(size_t n=0; n<min(lattice_result.size(), number_output); n++)
   	{
      	*os << "(" << ++count << ") ";

      	lattice_result[n].putReducedLatticeConstantsDegree(cdata.putBaseCenteredAxis(), cdata.putRhombohedralAxis(), length_axis, angle_axis);

      	*os << length_axis[0] << " " << length_axis[1] << " " << length_axis[2] << " "
      		<<  angle_axis[0] << " " <<  angle_axis[1] << " " <<  angle_axis[2] << " ("
      		<< put_bravais_type_name( lattice_result[n].enumBravaisType(), cdata.putBaseCenteredAxis() )
      	    << ", Volume=" << lattice_result[n].putLatticeVolume()
      		<< ", " << lattice_result[n].putFiguresOfMerit().putLabel_FigureOfMeritWolff()
      		<< "=" << lattice_result[n].putFiguresOfMerit().putFigureOfMeritWolff() << ")\n";
   	}
   	*os << endl;
}

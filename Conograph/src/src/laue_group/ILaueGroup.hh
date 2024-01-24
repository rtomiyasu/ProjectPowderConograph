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
#ifndef ILAUEGROUP_HH_
#define ILAUEGROUP_HH_

#include "../RietveldAnalysisTypes.hh"
#include "../utility_data_structure/index_set.hh"
#include "../utility_data_structure/VecDat3.hh"
#include "../point_group/enumPointGroup.hh"
#include "../symmetric_operation/MillerIndex.hh"
#include "../zerror_type/error_out.hh"

// Class of Laue Group.
class ILaueGroup
{
public:
	virtual string CrystalSystemName() const = 0;
	virtual string Name() const = 0;
	virtual Int4 LaueMultiplicity(const MillerIndex3&) const = 0;
	virtual void putLatticeConstantFlag(Mat_DP_constr&) const = 0;
	virtual bool checkLatticeConstant(VecDat3<Double>&, VecDat3<Double>&, Double thred = 1.0e-3) const = 0;
	inline ZErrorMessage checkLatticeConstantError(const VecDat3<Double>&, const VecDat3<Double>&, Double thred = 1.0e-3) const;
	virtual ~ILaueGroup(){};
};


inline string chToString(const VecDat3<Double>& length_axis, const VecDat3<Double>& angle_axis, const Int4& precision)
{
	return num2str(length_axis[0],precision) + ","
	 + num2str(length_axis[1],precision) + ","
	 + num2str(length_axis[2],precision) + ","
	 + num2str(angle_axis[0],precision) + ","
	 + num2str(angle_axis[1],precision) + ","
	 + num2str(angle_axis[2],precision);
}


inline ZErrorMessage ILaueGroup::checkLatticeConstantError(const VecDat3<Double>& length_axis,
		const VecDat3<Double>& angle_axis, Double thred) const
{
	VecDat3<Double> length_copy = length_axis;
	VecDat3<Double> angle_copy = angle_axis;

	if( !checkLatticeConstant(length_copy, angle_copy, thred ) )
	{
		return ZErrorMessage(ZErrorOutRange,
					"Lattice constants are out of range : "
					+chToString(length_axis, angle_axis, 6)
					+"\nProposed modification : "
					+chToString(length_copy, angle_copy, 6),
					__FILE__, __LINE__, __FUNCTION__);
	}
	return ZErrorMessage();
};

#endif /*ILAUEGROUP_H_*/

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
#ifndef enumBravaisType_HH_
#define enumBravaisType_HH_

#include "enumAxis.hh"
#include "../centring_type/enumCentringType.hh"
#include "../point_group/enumPointGroup.hh"

// Triclinic, Monoclinic, Monoclinic(C),
// Orthorhombic(P), Orthorhombic(C), Orthorhombic(I), Orthorhombic(F),
// Tetragonal(P), Tetragonal(I), Rhombohedral, Hexagonal.
// Cubic(P), Cubic(I), Cubic(F).
enum eBravaisType{ CSUndefined = -1, Triclinic = 0, Monoclinic_P = 1, Monoclinic_B = 2,
			Orthorhombic_P = 3, Orthorhombic_C = 4, Orthorhombic_I = 5, Orthorhombic_F = 6,
			Tetragonal_P = 7, Tetragonal_I = 8, Rhombohedral = 9, Hexagonal = 10,
			Cubic_P = 11, Cubic_I = 12, Cubic_F = 13 };

inline const Int4& put_number_of_bravais_types()
{
	static const Int4 NUM_LS = 14;
	return NUM_LS;
}


inline string put_bravais_type_name(const eBravaisType& i, const eABCaxis& axis)
{
	static const size_t NUM_LS = 14;
	static const string name[NUM_LS] = { "Triclinic", "Monoclinic(P)", "Monoclinic",
											"Orthorhombic(P)", "Orthorhombic(C)",
											"Orthorhombic(I)", "Orthorhombic(F)",
											"Tetragonal(P)", "Tetragonal(I)",
											"Rhombohedral", "Hexagonal",
											"Cubic(P)", "Cubic(I)", "Cubic(F)" };

	if( i == Monoclinic_B )
	{
		if( axis == A_Axis ) return name[(size_t)i] + "(B)";
		else if( axis == B_Axis ) return name[(size_t)i] + "(C)";
		else // if( axis == C_Axis )
			return name[(size_t)i] + "(A)";
	}
	else return name[(size_t)i];
}

#endif /*enumBravaisType_HH_*/

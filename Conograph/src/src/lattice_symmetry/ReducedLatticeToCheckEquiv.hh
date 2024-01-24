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
#ifndef LATTICESYMMETRY_HH_
#define LATTICESYMMETRY_HH_


#include "../lattice_symmetry/check_equiv.hh"


class ReducedLatticeToCheckEquiv
{
private:
	const Double m_resol;
	
protected:
	// Selling-reduced.
	const SymMat<Double> m_S_super;

	// Almost Selling-reduced matrix that is equivalent with m_S_super_obtuse.
	vector< SymMat<Double> > m_S_super_equiv;
	
	// Set m_S_super_equiv from m_S_super.
	void put_S_super_equiv(const Double& resol);

public:
	// S_super: Delaunay reduced.
	ReducedLatticeToCheckEquiv(const Double& resol,
									const SymMat<Double>& S_super);
	~ReducedLatticeToCheckEquiv();

	// Returns true if the argument S_super almost equivalent with m_S_super.
	bool equiv(const SymMat<Double>& S_super) const;
};

#endif /*LATTICESYMMETRY_HH_*/

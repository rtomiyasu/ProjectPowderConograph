/*
 * The MIT License

   BLDConograph (Bravais lattice determination module in Conograph)

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
#ifndef LATTICEFIGUREOFMERITTOCHECKSYMMETRY_HH_
#define LATTICEFIGUREOFMERITTOCHECKSYMMETRY_HH_

#include "LatticeFigureOfMerit.hh"

class OutputInfo;

// Class for outputting information about a lattice in index file.
class LatticeFigureOfMeritToCheckSymmetry
{
private:
	LatticeFigureOfMerit m_latfom;

	// m_S_red.first is Buerger-reduced.
	// m_S_red.second * m_S_red.first * Transpose(m_S_red.second) is Selling-reduced.
	SymMat43_Double m_S_red;

	bool checkIfLatticeIsMonoclinic(const ePointGroup& epg_new, const Double& resol,
								map< SymMat<Double>, NRMat<Int4> >& ans) const;
	bool checkIfLatticeIsOrthorhombic(const Double& resol,
								map< SymMat<Double>, NRMat<Int4> >& ans) const;
	bool checkIfLatticeIsTetragonal(const Double& resol,
								map< SymMat<Double>, NRMat<Int4> >& ans) const;
	bool checkIfLatticeIsHexagonal(const ePointGroup& epg_new, const Double& resol,
								map< SymMat<Double>, NRMat<Int4> >& ans) const;

public:
	LatticeFigureOfMeritToCheckSymmetry();
	LatticeFigureOfMeritToCheckSymmetry(const Double& rhs);	// Sets only m_latfom.m_determ_GramMat = rhs;
	LatticeFigureOfMeritToCheckSymmetry(const BravaisType& ebrat,
										const SymMat43_Double& S_red);
	virtual ~LatticeFigureOfMeritToCheckSymmetry(){};

	void setDistance(const SymMat<Double>& S_super0);
	inline bool operator<(const LatticeFigureOfMeritToCheckSymmetry& rhs) const { return m_latfom < rhs.m_latfom; };

	const LatticeFigureOfMerit& putLatticeFigureOfMerit() const { return m_latfom; };

	inline void setDeWolffFigureOfMerit(const Int4& num_ref_figure_of_merit, const vector<QData>& qdata) { m_latfom.setDeWolffFigureOfMerit(num_ref_figure_of_merit, -1, qdata); };

	// Returns true if the lattice has at least the symmetry of eps. 
	// On output, ans equals the equivalent lattice with symmetry of eps.
	bool checkLatticeSymmetry(const ePointGroup& epg_new, const Double& resol,
								map< SymMat<Double>, NRMat<Int4> >& ans) const;

	// Set-functions.
	void setLatticeConstants43(const BravaisType& brat, const SymMat43_Double& S);

	inline const SymMat43_Double& putInitialForm() const { return m_S_red; };
	inline const SymMat<Double> putInitialSellingReducedForm() const { return transform_sym_matrix(m_S_red.second, m_S_red.first); };

	void putLatticesOfHigherSymmetry(const ePointGroup& epg, const Double& resol,
									vector<LatticeFigureOfMeritToCheckSymmetry>& lattice_result) const;

	static bool cmpFOMdeWolff(const LatticeFigureOfMeritToCheckSymmetry& lhs, const LatticeFigureOfMeritToCheckSymmetry& rhs)
	{
		return LatticeFigureOfMerit::cmpFOMdeWolff( lhs.putLatticeFigureOfMerit(), rhs.putLatticeFigureOfMerit() );
	}
};

#endif /*LatticeFigureOfMeritToCheckSymmetry_HH_*/

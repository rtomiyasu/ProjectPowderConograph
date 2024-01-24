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
#ifndef LATTICEMETRICTENSOR_HH_
#define LATTICEMETRICTENSOR_HH_

#include "../utility_data_structure/nrutil_nr.hh"
#include "../utility_data_structure/SymMatNplus1N.hh"
#include "../utility_func/transform_sym_matrix.hh"
#include "../utility_lattice_reduction/matrix_NbyN.hh"

class HKL_Q;

// Class for outputting information about a lattice in index file.
class LatticeMetricTensor
{
	friend inline bool operator<(const LatticeMetricTensor& lhs, const LatticeMetricTensor& rhs);

private:
	// This matrix is a result of optimization, therefore, it may be not reduced.
	// m_S is Minkowski-reduced.
	// m_S.second * m_S.first * Transpose(m_S.second) is Selling-reduced.
	SymMatNplus1N_Double m_S;
	Double m_determ_S;
	
	Double m_figures_of_merit;
	Int4 m_num_lattice_found;
	
public:
	LatticeMetricTensor(const Int4& arg1, const Double& arg2);	// Sets only the size of m_S to arg1 and m_determ_S to arg2.
	// This assumes that S.first is Minkowski-reduced and S.second * S.first * Transpose(S.second) is Selling-reduced.
	LatticeMetricTensor(const SymMatNplus1N_Double& S);
	virtual ~LatticeMetricTensor(){}; 
	
	void setFigureOfMeritWu(const vector<Double>& qobs_tray,
							const vector< vector<HKL_Q>::const_iterator>& closest_qcal_tray);

	// Set-functions.
	// This method assumes that S.first is Minkowski-reduced and S.second * S.first * Transpose(S.second) is Selling-reduced.
	void setLatticeConstants(const SymMatNplus1N_Double& S);

	// Put-functions (Returns a value.)
	inline const SymMatNplus1N_Double& putReducedForm() const { return m_S; };
	inline const SymMat<Double>& putMinkowskiReducedForm() const { return m_S.first; };
	inline SymMat<Double> putSellingReducedForm() const { return transform_sym_matrix(m_S.second, m_S.first); };

	inline const Double& putDeterminantOfGramMatrix() const { return m_determ_S; };
	inline const Double& putFiguresOfMerit() const { return m_figures_of_merit; };

	// Returns unit-cell parameters with other centrings using m_S.
	static bool cmpFOMWu(const LatticeMetricTensor& lhs, const LatticeMetricTensor& rhs)
	{
		return lhs.m_figures_of_merit > rhs.m_figures_of_merit;
	}

	inline Int4& putNumberOfLatticesInNeighborhood() { return m_num_lattice_found; };
	inline const Int4& putNumberOfLatticesInNeighborhood() const { return m_num_lattice_found; };

};


inline bool operator<(const LatticeMetricTensor& lhs, const LatticeMetricTensor& rhs)
{ 
	return lhs.m_determ_S < rhs.m_determ_S;
}

#endif /*LATTICEMETRICTENSOR_HH_*/

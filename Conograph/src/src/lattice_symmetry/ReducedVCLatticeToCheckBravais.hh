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
#ifndef LATTICESYMMETRYPRIMTRICL_HH_
#define LATTICESYMMETRYPRIMTRICL_HH_

#include <set>
#include "../utility_data_structure/SymMat43_VCData.hh"
#include "../centring_type/enumCentringType.hh"
#include "../bravais_type/BravaisType.hh"

class FracMat;

class ReducedVCLatticeToCheckBravais
{
private:
	static const vector< pair< NRMat<Int4>, FracMat > > m_trans_mat_red_F;
	static const vector< pair< NRMat<Int4>, FracMat > > m_trans_mat_red_I;
	static const vector< vector< pair< NRMat<Int4>, FracMat > > > m_trans_mat_red_rhom;	// rho-axis(quick search/prudent search), hex-axis(quick search/prudent search).
	static const vector< vector< pair< NRMat<Int4>, FracMat > > > m_trans_mat_red_base;	// a-axis(quick search/prudent search), b-axis(quick search), c-axis(quick search/prudent search).

	const BravaisType m_monoclinic_b_type;
	const BravaisType m_rhombohedral_type;

	const SymMat<VCData> m_S_super_obtuse;

	// In the following member, equivalent lattices with symmetry of each Bravais-lattice are set.
	// m_S_super_obtuse = m_S_red_*.second * m_S_red_*.first * Transpose(m_S_red_*.second).
	map< SymMat<VCData>, NRMat<Int4> > m_S_red_base;
	map< SymMat<VCData>, NRMat<Int4> > m_S_red_body;
	map< SymMat<VCData>, NRMat<Int4> > m_S_red_face;
	map< SymMat<VCData>, NRMat<Int4> > m_S_red_rhom;

	
	// On input, inv_flag = false indicates that S_super_obtuse_equiv is Selling-reduced,
	// and inv_flag = true indicates that Inverse(S_super_obtuse_equiv) is Selling-reduced.
	// In the former case, on output, S_red_body are symmetric matrices having a body-centered and Minkowski-reduced inverse.
	// In the latter case, on output, S_red_IF are symmetric matrices having a face-centered and Minkowski-reduced inverse.
	static void put_S_Buerger_reduced_IF(
			const Double& cv2, const SymMat<VCData>& S_super_obtuse,
			map< SymMat<VCData>, NRMat<Int4> >& S_red_IF,
			const bool& inv_flag);

	static void put_S_Buerger_reduced_rhom(
			const BravaisType& rhombohedral_type,
			const bool& does_prudent_sym_search,
			const Double& cv2, const SymMat<VCData>& S_super_obtuse,
			map< SymMat<VCData>, NRMat<Int4> >& S_red_rhomhex);

	// On input, S_super_obtuse_equiv is Selling-reduced.
	// On output, S_red_base are symmetric matrices having a base-centered and Minkowski-reduced inverse.
	static void put_S_Buerger_reduced_base(
			const BravaisType& monoclinic_b_type,
			const bool& does_prudent_sym_search,
			const Double& cv2, const SymMat<VCData>& S_super_obtuse,
			map< SymMat<VCData>, NRMat<Int4> >& S_red_base);

public:
	ReducedVCLatticeToCheckBravais(const eABCaxis& axis1,
									const eRHaxis& axis2,
									const bool& does_prudent_sym_search,
									const Double& cv2, const SymMat<VCData>& S_super);
	~ReducedVCLatticeToCheckBravais();

	// Returns true if the lattice has at least the symmetry of eblat, .i.e., if the corresponding member m_S_red_* is not empty. 
	// On output, mat_trays equals the equivalent lattices with symmetry of eblat (m_S_red_*).
	const map< SymMat<VCData>, NRMat<Int4> >& checkCentringType(const BravaisType& brat) const;
};

#endif /*LATTICESYMMETRYPRIMTRICL_HH_*/

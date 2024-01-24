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
#ifdef DEBUG
#include <iostream>
#endif
#include "../utility_data_structure/index_set.hh"
#include "../utility_func/transform_sym_matrix.hh"
#include "../utility_lattice_reduction/put_Selling_reduced_lattice.hh"
#include "../utility_lattice_reduction/put_Buerger_reduced_lattice.hh"
#include "../lattice_symmetry/LatticeFigureOfMeritToCheckSymmetry.hh"
#include "../lattice_symmetry/ReducedLatticeToCheckBravais.hh"
#include "../zerror_type/error_out.hh"
#include "BravaisLattice.hh"

const bool BravaisLattice::m_DoesPrudentSymSearch = false;

BravaisLattice::BravaisLattice()
{
	m_num_ref_figure_of_merit = 20;
	m_resol = sqrt( numeric_limits<double>::epsilon() );
}


BravaisLattice::~BravaisLattice()
{
}


void BravaisLattice::putCentringTypes(const ReducedLatticeToCheckBravais& RLCB,
		const BravaisType& brat,
		vector<LatticeFigureOfMeritToCheckSymmetry>& lattice_result)
{
	lattice_result.clear();
	
	const map< SymMat<Double>, NRMat<Int4> >& S_red_tray = RLCB.checkCentringType(brat);
	if( S_red_tray.empty() ) return;

	// The lattice of RLCB has at least the symmetry given by eblat.
	SymMat<Double> S_super(4);
	NRMat<Int4> trans_mat(4,3);

	for(map< SymMat<Double>, NRMat<Int4> >::const_iterator it=S_red_tray.begin(); it!=S_red_tray.end(); it++)
	{
		S_super = transform_sym_matrix(it->second, it->first);
		trans_mat = identity_matrix<Int4>(4);

		// S_super = trans_mat * it->second * it->first * Transpose(trans_mat * it->second) is Delone reduced.
		if( !put_Selling_reduced_dim_less_than_4(S_super, trans_mat) )
		{
			assert( false );
		}
		moveSmallerDiagonalLeftUpper(S_super, trans_mat);

		lattice_result.push_back( LatticeFigureOfMeritToCheckSymmetry( brat, SymMat43_Double(it->first, mprod(trans_mat, it->second) ) ) );
	}
}



static SymMat43_Double putTransformMatrixFromSellingReducedToBuergerReduced(
		const SymMat<Double>& S_super_obtuse)
{
	NRMat<Int4> trans_mat(4,3);

	// S_red = trans_mat * S_super_obtuse * transpose(trans_mat).
	SymMat<Double> S_red(3);
	putBuergerReducedMatrix(S_super_obtuse, false, S_red, trans_mat);

	return SymMat43_Double( S_red, put_transform_matrix_row3to4( Inverse3( trans_mat ) ) );
}



void BravaisLattice::putLatticeCandidatesForEachBravaisTypes(const SymMat<Double>& S_super,
					const eABCaxis& abc_axis,
					const eRHaxis& rh_axis)
{
	try{

	for(Int4 i=1; i<NUM_LS; i++)
	{
		m_result[i].clear();
	}
	vector<LatticeFigureOfMeritToCheckSymmetry>& lattice_result_tri = m_result[(size_t)Triclinic];
	vector<LatticeFigureOfMeritToCheckSymmetry>& lattice_result_mono_P = m_result[(size_t)Monoclinic_P];
	vector<LatticeFigureOfMeritToCheckSymmetry>& lattice_result_mono_B = m_result[(size_t)Monoclinic_B];
	vector<LatticeFigureOfMeritToCheckSymmetry>& lattice_result_ortho_P = m_result[(size_t)Orthorhombic_P];
	vector<LatticeFigureOfMeritToCheckSymmetry>& lattice_result_ortho_B = m_result[(size_t)Orthorhombic_C];
	vector<LatticeFigureOfMeritToCheckSymmetry>& lattice_result_ortho_I = m_result[(size_t)Orthorhombic_I];
	vector<LatticeFigureOfMeritToCheckSymmetry>& lattice_result_ortho_F = m_result[(size_t)Orthorhombic_F];
	vector<LatticeFigureOfMeritToCheckSymmetry>& lattice_result_tetra_P = m_result[(size_t)Tetragonal_P];
	vector<LatticeFigureOfMeritToCheckSymmetry>& lattice_result_tetra_I = m_result[(size_t)Tetragonal_I];
	vector<LatticeFigureOfMeritToCheckSymmetry>& lattice_result_rhom = m_result[(size_t)Rhombohedral];
	vector<LatticeFigureOfMeritToCheckSymmetry>& lattice_result_hex = m_result[(size_t)Hexagonal];
	vector<LatticeFigureOfMeritToCheckSymmetry>& lattice_result_cubic_P = m_result[(size_t)Cubic_P];
	vector<LatticeFigureOfMeritToCheckSymmetry>& lattice_result_cubic_I = m_result[(size_t)Cubic_I];
	vector<LatticeFigureOfMeritToCheckSymmetry>& lattice_result_cubic_F = m_result[(size_t)Cubic_F];

	const SymMat43_Double S_red = putTransformMatrixFromSellingReducedToBuergerReduced(S_super);
	LatticeFigureOfMeritToCheckSymmetry latFOM(BravaisType( pair<eCentringType, ePointGroup>(Prim, Ci) ), S_red);
	lattice_result_tri.push_back( latFOM );

	// Calculate figures of merit as triclinic
	const ReducedLatticeToCheckBravais RLCB(abc_axis, rh_axis, m_DoesPrudentSymSearch, m_resol, S_super);

	vector<LatticeFigureOfMeritToCheckSymmetry> latFOM_tray;

	BravaisLattice::putCentringTypes(RLCB, BravaisType( put_monoclinic_b_type(abc_axis) ), latFOM_tray);
	lattice_result_mono_B.insert(lattice_result_mono_B.end(), latFOM_tray.begin(), latFOM_tray.end());

	BravaisLattice::putCentringTypes(RLCB, BravaisType( pair<eCentringType, ePointGroup>(Inner, D2h) ), latFOM_tray);
	lattice_result_ortho_I.insert(lattice_result_ortho_I.end(), latFOM_tray.begin(), latFOM_tray.end());

	BravaisLattice::putCentringTypes(RLCB, BravaisType( pair<eCentringType, ePointGroup>(Face, D2h) ), latFOM_tray);
	lattice_result_ortho_F.insert(lattice_result_ortho_F.end(), latFOM_tray.begin(), latFOM_tray.end());

	BravaisLattice::putCentringTypes(RLCB, BravaisType( put_rhombohedral_type(rh_axis) ), latFOM_tray);
	lattice_result_rhom.insert(lattice_result_rhom.end(), latFOM_tray.begin(), latFOM_tray.end());

	latFOM.putLatticesOfHigherSymmetry(put_monoclinic_p_type(abc_axis), m_resol, latFOM_tray);
	lattice_result_mono_P.insert(lattice_result_mono_P.end(), latFOM_tray.begin(), latFOM_tray.end());

	latFOM.putLatticesOfHigherSymmetry(D2h, m_resol, latFOM_tray);
	lattice_result_ortho_P.insert(lattice_result_ortho_P.end(), latFOM_tray.begin(), latFOM_tray.end());

	for(Int4 i=1; i<NUM_LS; i++)
	{
		const Int4 num_lattice = m_result[i].size();

		for(Int4 index=0; index<num_lattice; index++)
		{
			LatticeFigureOfMeritToCheckSymmetry& latFOM0 = m_result[i][index];

			if( eBravaisType(i) == Monoclinic_P )
			{
				latFOM0.putLatticesOfHigherSymmetry(D6h, m_resol, latFOM_tray);
				lattice_result_hex.insert(lattice_result_hex.end(), latFOM_tray.begin(), latFOM_tray.end());
			}
			else if( eBravaisType(i) == Monoclinic_B )
			{
				latFOM0.putLatticesOfHigherSymmetry(D2h, m_resol, latFOM_tray);
				lattice_result_ortho_B.insert(lattice_result_ortho_B.end(), latFOM_tray.begin(), latFOM_tray.end());
			}
			else if( eBravaisType(i) == Orthorhombic_P )
			{
				latFOM0.putLatticesOfHigherSymmetry(D4h_Z, m_resol, latFOM_tray);
				lattice_result_tetra_P.insert(lattice_result_tetra_P.end(), latFOM_tray.begin(), latFOM_tray.end());

				latFOM0.putLatticesOfHigherSymmetry(Oh, m_resol, latFOM_tray);
				lattice_result_cubic_P.insert(lattice_result_cubic_P.end(), latFOM_tray.begin(), latFOM_tray.end());
			}
			else if( eBravaisType(i) == Orthorhombic_I )
			{
				latFOM0.putLatticesOfHigherSymmetry(D4h_Z, m_resol, latFOM_tray);
				lattice_result_tetra_I.insert(lattice_result_tetra_I.end(), latFOM_tray.begin(), latFOM_tray.end());

				latFOM0.putLatticesOfHigherSymmetry(Oh, m_resol, latFOM_tray);
				lattice_result_cubic_I.insert(lattice_result_cubic_I.end(), latFOM_tray.begin(), latFOM_tray.end());
			}
			else if( eBravaisType(i) == Orthorhombic_F )
			{
				latFOM0.putLatticesOfHigherSymmetry(Oh, m_resol, latFOM_tray);
				lattice_result_cubic_F.insert(lattice_result_cubic_F.end(), latFOM_tray.begin(), latFOM_tray.end());
			}
		}
	}
    }
    catch(bad_alloc& ball){
    	throw nerror(ball, __FILE__, __LINE__, __FUNCTION__);
    }
}


void BravaisLattice::execute(const SymMat<Double>& S_super,
					const vector<QData>& set_of_qvalues, 	// set of q-values of S_super.
//					const Double& resol,
					const eABCaxis& abc_axis,
					const eRHaxis& rh_axis)
{
	putLatticeCandidatesForEachBravaisTypes(S_super, abc_axis, rh_axis);

	for(Int4 i=NUM_LS-1; i>=0; i--)
	{
		for(vector<LatticeFigureOfMeritToCheckSymmetry>::iterator it=m_result[i].begin(); it!=m_result[i].end(); it++)
		{
			it->setDeWolffFigureOfMerit(m_num_ref_figure_of_merit, set_of_qvalues);
		}
	}
}


const LatticeFigureOfMeritToCheckSymmetry& BravaisLattice::putBravaisLattice() const
{
	assert( m_result[0].size() == 1 );
	vector<LatticeFigureOfMeritToCheckSymmetry>::const_iterator ans = m_result[0].begin();
	for(Int4 i=NUM_LS-1; i>0; i--)
	{
		if( !m_result[i].empty() )
		{
			for(vector<LatticeFigureOfMeritToCheckSymmetry>::const_iterator it=m_result[i].begin(); it!=m_result[i].end(); it++)
			{
				if( LatticeFigureOfMeritToCheckSymmetry::cmpFOMdeWolff( *it, *ans ))
				{
					ans = it;
				}
			}
		}
	}
	return *ans;
}

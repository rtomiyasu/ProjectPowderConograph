/*
 * The MIT License

   LatticeWithSameQ (a module in Conograph to obtain all the lattices with the same computed lines)

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
#include "../ControlParam.hh"
#include "../bravais_type/BravaisLattice.hh"
#include "../utility_data_structure/index_set.hh"
#include "../utility_data_structure/range.hh"
#include "../utility_func/transform_sym_matrix.hh"
#include "../utility_lattice_reduction/matrix_NbyN.hh"
#include "../utility_lattice_reduction/put_Selling_reduced_lattice.hh"
#include "../zerror_type/error_out.hh"
#include "../zlog/zlog.hh"
#include "../lattice_symmetry/LatticeFigureOfMeritToCheckSymmetry.hh"
#include "../lattice_symmetry/ReducedLatticeToCheckEquiv.hh"
#include "../lattice_symmetry/gather_q_of_Ndim_lattice.hh"
#include "Choleskydcmp.hh"
#include "LatticeMetricTensor.hh"
#include "LatticeWithSameQ.hh"

LatticeWithSameQ::LatticeWithSameQ()
{
}


LatticeWithSameQ::~LatticeWithSameQ()
{
}


// Set the member variables.
void LatticeWithSameQ::setParam(const ControlParam& cont, const SymMat<Double>& S_super)
{
	m_set_of_qvalues.clear();
    m_result.clear();

    m_resol = cont.putResolutionTocheckComputedLines();
    m_dimension_of_lattice = S_super.size() - 1;
	assert( m_dimension_of_lattice > 0 );

	const Double max_Q = cont.putMaxQTocheckComputedLines();
	const Int4 Msq = m_scale_of_superlattice * m_scale_of_superlattice;

	gatherQcal(S_super, max_Q*Msq, m_set_of_qvalues);

	m_set_of_qvalues.erase(m_set_of_qvalues.begin(), upper_bound(m_set_of_qvalues.begin(), m_set_of_qvalues.end(), 0.0 ) );

	vector<Double>::iterator it_end = m_set_of_qvalues.end()-1;
	for(; it_end>=m_set_of_qvalues.begin()+100; it_end--)
	{
		if( *(it_end-1)*(1.0+m_resol) < *it_end*(1.0-m_resol) ) break;
	}
	m_set_of_qvalues.erase(it_end + 1, m_set_of_qvalues.end());

#ifdef DEBUG
ZLOG_INFO( "Outputting q-values...\nNo., qvalue\n" );
	stringstream strstream;
	for(size_t k=0; k<m_set_of_qvalues.size(); k++)
	{
		strstream << k + 1 << "  " << m_set_of_qvalues[k] << endl;
	}
ZLOG_INFO( strstream.str() + "\n" );
#endif
}

// qobs_tray: a sorted sequence <q_0, ..., q_{M-1}> of positive real numbers,
// resol: relative resolution of q-values,
// Msq: a natural number,
// S:=(s_{ij}) : N-by-N symmetric matrix,
// 0 <= m <= n < 4:  integers indicating the algorithm is exploring the options for S(m, n),
// assigned_qmin, assigned_qmax: range of Msq*S(n,n) or Msq*(S(m,m)+S(n,n)+S(m,n)*2).
// Ans: an array of N-by-N positive-definite Minkowski reduced symmetric matrices
//      such that whose set of representations contains qobs_tray and is contained in 1/scale_of_sublattice*qobs_tray
static void lattice_with_same_q(
		const vector<Double>& qobs_tray,
		const Int4& Msq,
		SymMat<Double>& S,
		const Int4& m,
		const Int4& n,
		const vector<Double>::const_iterator& assigned_qmin,
		const vector<Double>::const_iterator& assigned_qmax,
		const vector<Double>::const_iterator& first_qobs_notin_qcal,
		const Double& resol,
		vector< SymMat<Double> >& ans,
		pair<bool, Double>& evaluation_ans)
{
	assert( m <= n );
	assert( m < S.size() );
	assert( n < S.size() );

	vector<Double>::const_iterator new_assigned_qmin, new_assigned_qmax, itM2;
	vector< vector<HKL_Q>::const_iterator > closest_qcal_tray;
	for(vector<Double>::const_iterator it=assigned_qmin; it<assigned_qmax; it++)
	{
		if ( m == n )
		{
			S(n, n) = *it / Msq;
		}
		else
		{
			S(m, n) = ( *it / Msq - S(m, m) - S(n, n) ) * 0.5;
		}

		if ( m <= 0 )
		{
			if ( n + 1 >= S.size() )
			{
				ans.push_back(S);
			}
			else
			{
				SymMat<Double> S_cp(n+1);
				for(Int4 i=0; i<n+1; i++)
				{
					for(Int4 j=0; j<=i; j++)
					{
						S_cp(i,j) = S(i,j);
					}
				}

				SymMat<Double> S_super(n+2);
				if( !put_Selling_reduced_matrix(S_cp, S_super) ) continue;

				vector<HKL_Q> qcal_tray;
				Double max_Q = min(S(n,n)*8.0, *(qobs_tray.rbegin())*(1.0+resol));
				gatherQcal(S_super, max_Q, qcal_tray);
				// qcal_tray*Msq is included in qobs_tray?
				if( !associateQcalWithQobs( upper_bound(qcal_tray.begin(), qcal_tray.end(), HKL_Q(NRVec<Int4>(), 0.0)),
											upper_bound(qcal_tray.begin(), qcal_tray.end(), HKL_Q(NRVec<Int4>(), *(qobs_tray.rbegin())/((1.0+resol)*Msq))),
											Msq, qobs_tray, resol) ) continue;

				// first entry of qobs_tray not contained in qcal_tray.
				vector<Double>::const_iterator first_qobs_notin_qcal2 = first_qobs_notin_qcal;
				do
				{
					itM2 = associateQobsWithQcal(first_qobs_notin_qcal2, qobs_tray.end(), qcal_tray, resol, closest_qcal_tray);

					if( itM2 == qobs_tray.end() )
					{
						evaluation_ans.first = false;
						break;
					}
					if( *itM2 < max_Q )
					{
						break;
					}
					if( max_Q >= *(qobs_tray.rbegin())*(1.0+resol) )
					{
						evaluation_ans.first = false;
						break;
					}

					max_Q = min(max_Q*2.0, *(qobs_tray.rbegin())*(1.0+resol));
					first_qobs_notin_qcal2 = itM2;
					gatherQcal(S_super, max_Q, qcal_tray);
				} while ( true );

				new_assigned_qmin = closest_data( qobs_tray.begin(), qobs_tray.end(), S(n,n)*Msq );
				new_assigned_qmax = (itM2== qobs_tray.end()?
										itM2:closest_data( qobs_tray.begin(), qobs_tray.end(), *itM2*Msq )+1);
				assert( new_assigned_qmin < new_assigned_qmax );

				evaluation_ans.second = max(evaluation_ans.second, *itM2*Msq);

				lattice_with_same_q(qobs_tray, Msq, S, n+1, n+1,
						new_assigned_qmin, new_assigned_qmax, itM2, resol, ans, evaluation_ans);
			}
		}
		else
		{
			Double Sm_1n_min, Sm_1n_max;
			Sm_1n_min = -S(m-1, m-1);
			if( m == n )
			{
				Sm_1n_max = 0;
			}
			else // 1 <= m < n
			{
				Sm_1n_max = S(m-1, m-1);

				// Case of m < n.
				// -2( S(m-1,n) + S(m,n) ) <= S(m-1,m-1) + S(m,m) + S(m-1,m)*2.
				Sm_1n_min = max( Sm_1n_min, - S(m-1,m-1) - S(m,m) - (S(m-1,m) + S(m,n))*2.0 );

				if( m < n - 1 ) // (m-1, n) = (0, 3)
				{
					// 2 | S(0,3) + S(1,3) | <= S(0,0) + S(1,1) + 2 S(0,1),
					const Double q01 = S(0,0) + S(1,1) + S(0,1)*2.0;
					Sm_1n_min = max( Sm_1n_min, -q01 - S(1,3)*2.0 );
					Sm_1n_max = min( Sm_1n_max,  q01 - S(1,3)*2.0 );

					// 2 | S(0,2) + S(0,3) | <= S(0,0) + S(2,2) + 2 S(2,3),
					const Double q02_23 = S(0,0) + S(2,2) + S(2,3)*2.0;
					Sm_1n_min = max( Sm_1n_min, -q02_23 - S(0,2)*2.0 );
					Sm_1n_max = min( Sm_1n_max,  q02_23 - S(0,2)*2.0 );

					// 2 | S(0,3) + S(1,3) + S(2,3) | <= S(0,0) + S(1,1) + S(2,2) + 2 S(0,1) + 2 S(0,2) + 2 S(1,2),
					const Double q012 = S(0,0) + S(1,1) + S(2,2) + (S(0,1) + S(0,2) + S(1,2))*2.0;
					Sm_1n_min = max( Sm_1n_min, -q012 - (S(1,3) + S(2,3))*2.0 );
					Sm_1n_max = min( Sm_1n_max,  q012 - (S(1,3) + S(2,3))*2.0 );

					// 2 | S(0,3) - S(1,3) - S(2,3) | <= S(0,0) + S(1,1) + S(2,2) - 2 S(0,1) - 2 S(0,2) + 2 S(1,2),
					const Double q_012 = S(0,0) + S(1,1) + S(2,2) - (S(0,1) + S(0,2) - S(1,2))*2.0;
					Sm_1n_min = max( Sm_1n_min, -q_012 + (S(1,3) + S(2,3))*2.0 );
					Sm_1n_max = min( Sm_1n_max,  q_012 + (S(1,3) + S(2,3))*2.0 );

					// 2 ( S(0,3) + S(1,3) - S(2,3) ) <= S(0,0) + S(1,1) + S(2,2) + 2 S(0,1) - 2 S(0,2) - 2 S(1,2).
					Sm_1n_max = min( Sm_1n_max,  S(0,0) + S(1,1) + S(2,2) + (S(0,1) - S(0,2) - S(1,2) - S(1,3) + S(2,3))*2.0 );
				}
			}

			new_assigned_qmin = lower_bound( qobs_tray.begin(), qobs_tray.end(), (S(m-1,m-1)+S(n,n)+Sm_1n_min)*Msq*(1.0 - resol) );
			new_assigned_qmax = upper_bound( new_assigned_qmin, qobs_tray.end(), (S(m-1,m-1)+S(n,n)+Sm_1n_max)*Msq*(1.0 + resol) );

			if( (S(m-1,m-1)+S(n,n)+Sm_1n_max)*Msq > *(qobs_tray.rbegin()) )
			{
				evaluation_ans.first = false;
			}

			lattice_with_same_q(qobs_tray, Msq, S, m-1, n,
				new_assigned_qmin, new_assigned_qmax, first_qobs_notin_qcal, resol, ans, evaluation_ans);
		}
	}
}

// On output,
// evaluation_ans.first: all the solusions were gained?
// evaluation_ans.second: max of q-values used to determine the diagonal entries of solutions.
// The second argument arg is used in the member method setNumberOfNeighbors.
void LatticeWithSameQ::execute(pair<bool, Double>& evaluation_ans, const SymMat<Double>& arg)
{
	m_result.clear();
	if( m_set_of_qvalues.empty() ) return;

	const Double Msq = m_scale_of_superlattice*m_scale_of_superlattice;
	const vector<Double>::const_iterator assigned_qmax
		= closest_data( m_set_of_qvalues.begin(), m_set_of_qvalues.end(), *(m_set_of_qvalues.begin())*Msq ) + 1;

	SymMat<Double> S(m_dimension_of_lattice);
	vector< SymMat<Double> > ans;

	evaluation_ans.first = true;
	evaluation_ans.second = 0.0;
	lattice_with_same_q(m_set_of_qvalues, Msq, S, 0, 0,
						m_set_of_qvalues.begin(), assigned_qmax, m_set_of_qvalues.begin(), m_resol, ans, evaluation_ans);

#ifdef DEBUG
	ZLOG_INFO( "Number of generated solutions: " + num2str(ans.size()) + "\n\n" );
#endif

	SymMat<Double> S_super(m_dimension_of_lattice+1);
	NRMat<Int4> TransMat(m_dimension_of_lattice+1, m_dimension_of_lattice);
	vector<HKL_Q> qcal_tray;
	vector< vector<HKL_Q>::const_iterator > closest_qcal_tray;

	for(vector< SymMat<Double> >::iterator it=ans.begin(); it<ans.end(); it++)
	{
		if( !put_Selling_reduced_dim_less_than_4(*it, S_super, TransMat) ) continue;

		gatherQcal(S_super, *(m_set_of_qvalues.rbegin())*(1.0+m_resol), qcal_tray);

		if( associateQcalWithQobs(upper_bound(qcal_tray.begin(), qcal_tray.end(), HKL_Q(NRVec<Int4>(), 0.0)),
									upper_bound(qcal_tray.begin(), qcal_tray.end(), HKL_Q(NRVec<Int4>(), *(m_set_of_qvalues.rbegin())/((1.0+m_resol)*Msq))),
									Msq, m_set_of_qvalues, m_resol)
			&& m_set_of_qvalues.end() == associateQobsWithQcal(m_set_of_qvalues.begin(), m_set_of_qvalues.end(), qcal_tray, m_resol, closest_qcal_tray) )
		{
			moveSmallerDiagonalLeftUpper(S_super, TransMat);

			LatticeMetricTensor LMT( SymMatNplus1N_Double(*it, TransMat) );
			LMT.setFigureOfMeritWu(m_set_of_qvalues, closest_qcal_tray);

			m_result.push_back( LMT );
		}
	}

#ifdef DEBUG
	ZLOG_INFO( "Number of solutions satisfying the required properties: " + num2str(m_result.size()) + "\n" );
	ZLOG_INFO( "Removing almost equivalent solutions...\n" );
#endif

	this->setNumberOfNeighbors(arg);

}

void LatticeWithSameQ::setNumberOfNeighbors(const SymMat<Double>& S_super0)
{
	stable_sort( m_result.begin(), m_result.end() ); // Sort by the unit-cell volume.
	for(vector<LatticeMetricTensor>::iterator it=m_result.begin(); it<m_result.end(); it++)
	{
		it->putNumberOfLatticesInNeighborhood() = 0;
	}

	const Double coef_lower = 1.0 - m_resol*3.0;
	const Double coef_upper = 1.0 + m_resol*3.0;
	if( S_super0.size() == m_dimension_of_lattice + 1 )
	{
		const Double detS0 = Determinant( put_sym_matrix_sizeNplus1toN(S_super0) );
		const vector<LatticeMetricTensor>::iterator it_begin = lower_bound( m_result.begin(), m_result.end(), LatticeMetricTensor(m_dimension_of_lattice, detS0*coef_lower) );
		const vector<LatticeMetricTensor>::iterator it_end = upper_bound( it_begin, m_result.end(), LatticeMetricTensor(m_dimension_of_lattice, detS0*coef_upper ) );

		const ReducedLatticeToCheckEquiv RLCS0(m_resol, S_super0);
		for(vector<LatticeMetricTensor>::iterator it=it_begin; it<it_end; it++)
		{
			if( RLCS0.equiv( it->putSellingReducedForm() ) )
			{
				it->putNumberOfLatticesInNeighborhood() = -1;
			}
		}
	}

	const Int4 num_lattice = m_result.size();
	for(Int4 index=0; index<num_lattice; index++)
	{
		LatticeMetricTensor& latFOM0 = m_result[index];
		if( latFOM0.putNumberOfLatticesInNeighborhood() < 0 ) continue;

		const Double& detS = latFOM0.putDeterminantOfGramMatrix();
		const Int4 ibegin = distance( m_result.begin(), lower_bound( m_result.begin(), m_result.end(), LatticeMetricTensor(m_dimension_of_lattice, detS*coef_lower) ) );
		const Int4 iend = distance( m_result.begin(), upper_bound( m_result.begin() + ibegin, m_result.end(), LatticeMetricTensor(m_dimension_of_lattice, detS*coef_upper ) ) );

		Int4 count=0;
		const ReducedLatticeToCheckEquiv RLCS(m_resol, latFOM0.putSellingReducedForm());
		for(Int4 index2=ibegin; index2<iend; index2++)
		{
			if( index2 == index ) continue;

			LatticeMetricTensor& latFOM2 = m_result[index2];

			// m_result_tri[index2] equals trans_mat * RLCB.m_S_super_obtuse * Transpose(trans_mat)
			if( RLCS.equiv( latFOM2.putSellingReducedForm() ) )
			{
				// Compare the figures of merit.
				if( LatticeMetricTensor::cmpFOMWu( latFOM2, latFOM0 ) )
				{
					if( latFOM2.putNumberOfLatticesInNeighborhood() >= 0 )
					{
						latFOM2.putNumberOfLatticesInNeighborhood() += 1;
						count = -1;
						break;
					}
				}
				else
				{
					count++;
					latFOM2.putNumberOfLatticesInNeighborhood() = -1;
				}
			}
		}

		latFOM0.putNumberOfLatticesInNeighborhood() = count;
	}

	Int4 index = 0;
	for(vector<LatticeMetricTensor>::const_iterator it=m_result.begin(); it<m_result.end(); it++)
   	{
		if( it->putNumberOfLatticesInNeighborhood() >= 0 ) index++;
   	}
}


void LatticeWithSameQ::putLatticesWithSameQ(const ControlParam& cdata,
		const vector<QData>& set_of_qvalues,
		vector<LatticeFigureOfMerit>& tray) const
{
	tray.clear();
	for(vector<LatticeMetricTensor>::const_iterator it=m_result.begin(); it<m_result.end(); it++)
   	{
		if( it->putNumberOfLatticesInNeighborhood() < 0 ) continue;

		BravaisLattice brav_lat;
    	brav_lat.setParam(cdata.putNumberOfReflectionsForFigureOfMerit(), m_resol);
      	brav_lat.execute(it->putSellingReducedForm(), set_of_qvalues, // m_resol,
      						cdata.putBaseCenteredAxis(), cdata.putRhombohedralAxis() );

      	tray.push_back( brav_lat.putBravaisLattice().putLatticeFigureOfMerit() );
   	}
}

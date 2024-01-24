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
#ifdef DEBUG
#include <iostream>
#endif
#include <limits>
#include "../utility_lattice_reduction/matrix_NbyN.hh"
#include "../utility_lattice_reduction/put_Selling_reduced_lattice.hh"
#include "check_equiv.hh"
#include "ReducedLatticeToCheckEquiv.hh"


ReducedLatticeToCheckEquiv::ReducedLatticeToCheckEquiv(const Double& resol,
		const SymMat<Double>& S_super)
	: 	m_resol(resol),
		m_S_super( S_super )
{
	assert( 1 < S_super.size() && S_super.size() <= 5 );
#ifdef DEBUG
	const Int4 ISIZE = S_super.size();
	Int4 count = 0;
	for(Int4 i=0; i<ISIZE; i++)
	{
		for(Int4 j=0; j<i; j++)
		{
			if( S_super(i,j) > 0.0 ) count++;
		}
	}
	if( ISIZE == 5 )
	{
		assert( count < 3 );
	}
	else assert( count <= 0 );
#endif
	this->put_S_super_equiv(resol);
}


ReducedLatticeToCheckEquiv::~ReducedLatticeToCheckEquiv()
{
}


bool ReducedLatticeToCheckEquiv::equiv(const SymMat<Double>& S_super) const
{
	assert( S_super.size() == m_S_super.size() );

	for(vector< SymMat<Double> >::const_iterator it=m_S_super_equiv.begin(); it!=m_S_super_equiv.end(); it++)
	{
		if( check_equiv_s(*it, S_super, m_resol) )
		{
			return true;
		}
	}
	return false;
}


static bool operator<(const SymMat<Double>& lhs, const SymMat<Double>& rhs)
{
	static const Double EPS_1 = 1.0+sqrt( numeric_limits<double>::epsilon() );
	assert( lhs.size() == rhs.size() );

	const Int4 ISIZE = lhs.size();
	for(Int4 i=0; i<ISIZE; i++)
	{
		for(Int4 j=0; j<i; j++)
		{
			const Double lhs_ij = lhs(i,i)+lhs(j,j)+lhs(i,j)*2.0;
			const Double rhs_ij = rhs(i,i)+rhs(j,j)+rhs(i,j)*2.0;

			if( lhs_ij*EPS_1 < rhs_ij ) return true;
			if( rhs_ij*EPS_1 < lhs_ij ) return false;
		}
	}

	return false;
}



inline Double put_diagonal_sum(const SymMat<Double>& S)
{
	const Int4 ISIZE = S.size();
	Double sum = 0.0;
	for(Int4 i=0; i<ISIZE; i++)
	{
		sum += S(i,i);
	}
	return sum;
}

static void insert_super_equiv_dim_4(const SymMat<Double>& S_super,
		const Double& resol, const Double& Max_totalQ,
		vector< SymMat<Double> >& equivalent_tray)
{
	const Int4 ISIZE = S_super.size();
	assert( ISIZE == 4 );

	const Double totalQ = put_diagonal_sum(S_super);
	const Double coef = 2.0;
#ifdef DEBUG
	if( totalQ <= 0.0 )
	{
		stringstream strstream;
		for(Int4 k=0; k<S_super.size(); k++)
		{
			for(Int4 k2=0; k2<=k; k2++)
			{
				strstream << S_super(k,k2) << "  ";
			}
			strstream << endl;
		}
		cout << strstream.str() + "\n";
	}
#endif

	assert( totalQ > 0.0 );

	Int4 k,l;
	vector< SymMat<Double> >::iterator it_begin;
	for(Int4 i=0; i<ISIZE-1; i++)
	{
		for(Int4 j=i+1; j<ISIZE; j++)
		{
			const Double& Sij = S_super(i,j);

			if( 0.0 < Sij ) continue;
			if( !equiv_zero(S_super, i,j, resol) ) continue;
			if( Max_totalQ + Sij * coef < totalQ ) continue;

			const NRMat<Int4>& mat = put_reduction_matrix(ISIZE, i,j);
			SymMat<Double> new_S_super( transform_sym_matrix(mat, S_super) );

			put_complement_set4(i,j,k,l);
			// Ki, Kj, Kk, Kl -> -Ki, Kj, Ki+Kk, Ki+Kl => new_Skl = Skl - Sij,
			// sum_{i=0}^4 new_Sii = -2*Sij + sum_{i=0}^4 Sii.
			if( 0.0 < new_S_super(k,l) && !equiv_zero(new_S_super, k,l, resol) ) continue;

			moveSmallerDiagonalLeftUpper(new_S_super);
			it_begin = lower_bound(equivalent_tray.begin(), equivalent_tray.end(), new_S_super);
			if( it_begin == equivalent_tray.end() || new_S_super < *it_begin )
			{
				assert( it_begin == equivalent_tray.begin() || *(it_begin-1) < new_S_super );
				equivalent_tray.insert( it_begin, new_S_super );
				insert_super_equiv_dim_4( new_S_super, resol, Max_totalQ, equivalent_tray );
			}
		}
	}
};



static void set_super_equiv(const SymMat<Double>& S_super, const Double& resol,
		vector< SymMat<Double> >& equivalent_tray)
{
	assert( S_super.size() == 4 );
	equivalent_tray.clear();

	equivalent_tray.push_back( S_super );

	const Double Max_totalQ = put_diagonal_sum(S_super) * (1.0 + resol);
	insert_super_equiv_dim_4( S_super, resol, Max_totalQ, equivalent_tray );
};



inline void permute_row_column(const Int4& i, const Int4& j, SymMat<Double>& rhs)
{
	const NRMat<Int4> tmat = put_permutation_matrix(rhs.size(), i, j);
	rhs = transform_sym_matrix( tmat, rhs );
};


// On input, for any stage < i < j < S_super.size(), S_super(i, i) <= S_super(j, j)
static void insert_super_permuted(const Int4& stage,
		const SymMat<Double>& S_super,
		const Double& resol,
		vector< SymMat<Double> >& equivalent_tray)
{
	const Int4 ISIZE = S_super.size();
	if( stage + 1 >= ISIZE )
	{
		const vector< SymMat<Double> >::iterator it_begin = lower_bound(equivalent_tray.begin(), equivalent_tray.end(), S_super);
		if( it_begin == equivalent_tray.end() || S_super < *it_begin )
		{
			assert( it_begin == equivalent_tray.begin() || *(it_begin-1) < S_super );
			equivalent_tray.insert( it_begin, S_super );
		}
		return;
	}

	insert_super_permuted(stage+1, S_super, resol, equivalent_tray);

	SymMat<Double> S_super_permuted = S_super;

	for(Int4 i=stage+1; i<ISIZE; i++)
	{
		if( !equiv_resol(S_super(stage, stage), S_super(i,i), resol ) ) break;
		permute_row_column( stage, i, S_super_permuted );
		insert_super_permuted(stage+1, S_super_permuted, resol, equivalent_tray);
	}
};


void set_super_permuted(const SymMat<Double>& S_super,
		const Double& resol,
		vector< SymMat<Double> >& equivalent_tray)
{
	equivalent_tray.clear();
	insert_super_permuted( 0, S_super, resol, equivalent_tray );
};



void ReducedLatticeToCheckEquiv::put_S_super_equiv(const Double& resol)
{
	m_S_super_equiv.clear();
	if( m_S_super.size() <= 2 )
	{
		m_S_super_equiv.push_back( m_S_super );
		return;
	}

	vector< SymMat<Double> > super_equiv_tray;
	set_super_equiv( m_S_super, resol, super_equiv_tray );

	vector< SymMat<Double> > permute_equiv_tray;
	vector< SymMat<Double> >::iterator it_begin;
	for(vector< SymMat<Double> >::iterator it=super_equiv_tray.begin(); it!=super_equiv_tray.end(); it++)
	{
		set_super_permuted( *it, resol, permute_equiv_tray );

		for(vector< SymMat<Double> >::const_iterator it3=permute_equiv_tray.begin(); it3!=permute_equiv_tray.end(); it3++)
		{
			it_begin = lower_bound(m_S_super_equiv.begin(), m_S_super_equiv.end(), *it3);
			if( it_begin == m_S_super_equiv.end() || *it3 < *it_begin )
			{
				assert( it_begin == m_S_super_equiv.begin() || *(it_begin-1) < *it3 );
				m_S_super_equiv.insert(it_begin, *it3);
			}
		}
	}
}

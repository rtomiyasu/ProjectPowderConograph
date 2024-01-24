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
#include <limits>
#include <algorithm>
#include "gather_q_of_Ndim_lattice.hh"
#include "../utility_func/zmath.hh"


// S_super is a Selling-reduced positive definite symmetric matrix.
// - sum_{i != j} (ci-cj)^2*Sij <= maxQ
static void set_max_range(const SymMat<Double>& S_super, const Double& maxQ,
		SymMat<Int4>& max_range)
{
	const Int4 ISIZE = max_range.size();
	assert( S_super.size() == ISIZE + 1 );
	assert( 1 <= ISIZE && ISIZE <= 4 );

	static const Int4 MAX_INT = numeric_limits<Int4>::max();
	static const Int4 SQUARE_MAX_INT = ifloor( sqrt( (Double)MAX_INT ) );
	const Double min_dp = - maxQ / Double(MAX_INT);
	max_range = SQUARE_MAX_INT;

	// From (ci-cj)^2 <= (ci-ck)^2 + (cj-ck)^2 + (ci-cl)^2 + (cj-cl)^2,
	// Sij > 0
	// => -Sij*(ci-cj)^2 >= -Sij * { (ci-ck)^2 + (cj-ck)^2 + (ci-cl)^2 + (cj-cl)^2 }.
	// Sij > 0, Skl > 0 for distinct 1 <= i, j, k, l, m <= 5
	// => -Sij*(ci-cj)^2 - Skl*(ck-cl)^2 >= -(Sij + Skl) * { (ci-ck)^2 + (cj-ck)^2 + (ci-cl)^2 + (cj-cl)^2 },
	//    -Sij*(ci-cj)^2 - Skl*(ck-cl)^2 >= -Sij * { (ci-cm)^2 + (cj-cm)^2 + (ci-cl)^2 } - Skl * { (cm-ck)^2 + (cj-ck)^2 + (cm-cl)^2 } - (Sij + Skl) * (cj-cl)^2
	for(Int4 n=0; n<ISIZE; n++)
	{
		for(Int4 n2=n+1; n2<ISIZE; n2++)
		{
			if( S_super(n, n2) < min_dp ) max_range(n, n2) = ifloor( sqrt( - maxQ / S_super(n, n2) ) );
		}
		if( S_super(n, ISIZE) < min_dp ) max_range(n, n) = ifloor( sqrt( - maxQ / S_super(n, ISIZE) ) );
	}
}


static void arrange_max_range(SymMat<Int4>& max_range)
{
	const Int4 ISIZE = max_range.size();

	Int4 num;
	bool flag = true;
	while( flag )
	{
		flag = false;
		
		for(Int4 i=0; i<ISIZE; i++)
		{
			for(Int4 i2=0; i2<i; i2++)
			{
				num = max_range(i,i)+max_range(i2,i2);
				for(Int4 k=0; k<ISIZE; k++)
				{
					// |ci - ci2| <= |ci - ck| + |ci2 - ck|
					num = min( max_range(i,k)+max_range(i2,k), num);
				}
				if( num < max_range(i,i2) )
				{
					max_range(i,i2) = num;
					flag = true;
				}
			}

			num = max_range(i,i);
			for(Int4 k=0; k<ISIZE; k++)
			{
				// |ci| <= |ci - ck| + |ck|
				num = min( max_range(i,k)+max_range(k,k), num );
			}
			if( num < max_range(i,i) )
			{
				max_range(i,i) = num;
				flag = true;
			}
			
		}
	}
}


inline Double norm(const NRVec<Int4>& vec_Zn, const SymMat<Double>& S_super)
{
	const Int4 ISIZE = vec_Zn.size();
	assert( ISIZE + 1 == S_super.size() );

	Double ans = 0.0;
	for(Int4 i=0; i<ISIZE; i++)
	{
		for(Int4 j=0; j<i; j++)
		{
			ans += S_super(i,j)*(vec_Zn[i]*vec_Zn[j]*2);
		}
		ans += S_super(i,i)*(vec_Zn[i]*vec_Zn[i]);
	}
	return ans;
}

static void set_candidate_Q(const Int4& index,
		const SymMat<Double>& S_super,
		const Double& maxQ,
		const SymMat<Int4>& max_range, 
		NRVec<Int4>& vec_Zn, vector<HKL_Q>& qcal_tray)
{
	const Int4 ISIZE = vec_Zn.size();
	assert( ISIZE + 1 == S_super.size() );

	if( index >= ISIZE )	// coef is determined.
	{
		const Double Q = norm( vec_Zn, S_super );
		if( maxQ >= Q )
		{
			qcal_tray.push_back(HKL_Q(vec_Zn, Q));
		}
		return;
	}
	
	// vec_Zn[index] <= max_coef
	// vec_Zn[index] >= min_coef
	Int4 max_coef = max_range(index,index);
	Int4 min_coef = -max_range(index,index);
	for(Int4 k=0; k<index; k++)
	{
		// vec_Zn[index] - vec_Zn[k] <= max_range(k,index).
		max_coef = min( max_coef, max_range(k,index) + vec_Zn[k] );

		// -max_range(k,index) <= vec_Zn[index] - vec_Zn[k].
		min_coef = max( min_coef, -max_range(k,index) + vec_Zn[k] );
	}
	
	// First non-zero entry should be positive.
	bool non_zero_entry = false;
	for(Int4 i=0; i<index; i++)
	{
		if( vec_Zn[i] != 0 )
		{
			non_zero_entry = true;
			break;
		}
	}
	if( !non_zero_entry ) min_coef = 0;

	for(Int4 ic=max_coef; ic>=min_coef; ic--)
	{
		vec_Zn[index] = ic;
		set_candidate_Q(index+1, S_super, maxQ, max_range, vec_Zn, qcal_tray);
	}
}


void gatherQcal(const SymMat<Double>& S_super, 
		const Double& maxQ,
		vector<HKL_Q>& qcal_tray
	)
{
	qcal_tray.clear();

	const Int4 ISIZE = S_super.size() - 1;
	assert( 1 <= ISIZE && ISIZE <= 4 );

	SymMat<Int4> max_range(ISIZE);
	set_max_range(S_super, maxQ, max_range);
	arrange_max_range(max_range);

	NRVec<Int4> vec_ZN(ISIZE);
	set_candidate_Q(0, S_super, maxQ, max_range, vec_ZN, qcal_tray);

	sort( qcal_tray.begin(), qcal_tray.end() );
}


inline VecDat3<Int4> product_hkl(const VecDat3<Int4>& lhs, const NRMat<Int4>& rhs)
{
	assert( rhs.nrows() >= 3 && rhs.ncols() == 3 );

	VecDat3<Int4> ans;
	ans[0] = lhs[0]*rhs[0][0] + lhs[1]*rhs[1][0] + lhs[2]*rhs[2][0];
	ans[1] = lhs[0]*rhs[0][1] + lhs[1]*rhs[1][1] + lhs[2]*rhs[2][1];
	ans[2] = lhs[0]*rhs[0][2] + lhs[1]*rhs[1][2] + lhs[2]*rhs[2][2];
	return ans;
}


// On input, S_super = TransMat * S * transpose(TransMat).
void gatherQcal(const SymMat<Double>& S_super,
		const Double& maxQ,
		const NRMat<Int4>& transform_hkl,
		vector<HKL_Q>& qcal_tray
	)
{
	gatherQcal(S_super, maxQ, qcal_tray);

	for(vector<HKL_Q>::iterator it=qcal_tray.begin(); it<qcal_tray.end(); it++)
	{
		it->setHKL( product_hkl(it->HKL(), transform_hkl) );
	}
}


bool associateQcalWithQobs(
		const vector<HKL_Q>::const_iterator& it_begin,
	    const vector<HKL_Q>::const_iterator& it_end,
		const Int4& scale_of_qcal,
		const vector<Double>& qobs_tray,
		const Double& resol)
{
	vector<Double>::const_iterator it_begin2, it_end2;
	for(vector<HKL_Q>::const_iterator it=it_begin; it<it_end; it++)
	{
		it_begin2 = lower_bound( qobs_tray.begin(), qobs_tray.end(), it->Q()*scale_of_qcal*(1.0 - resol) );
		it_end2 = upper_bound( it_begin2, qobs_tray.end(), it->Q()*scale_of_qcal*(1.0 + resol) );
		if( it_begin2 >= it_end2 ) return false;
	}
	return true;
}


void associateQobsWithQcal(
		vector<Double>::const_iterator it_begin, vector<Double>::const_iterator it_end,
		const vector<HKL_Q>& qcal_tray,
		vector< vector<HKL_Q>::const_iterator >& closest_qcal_tray)
{
	closest_qcal_tray.clear();

	for(vector<Double>::const_iterator it=it_begin; it<it_end; it++)
	{
		closest_qcal_tray.push_back( closest_data(qcal_tray.begin(), qcal_tray.end(), *it) );
	}
}


vector<Double>::const_iterator associateQobsWithQcal(
			const vector<Double>::const_iterator& it_begin,
		    const vector<Double>::const_iterator& it_end,
		    const vector<HKL_Q>& qcal_tray, const Double& resol,
		    vector< vector<HKL_Q>::const_iterator >& closest_qcal_tray)
{
	closest_qcal_tray.clear();

	vector<HKL_Q>::const_iterator it_begin2, it_end2;
	for(vector<Double>::const_iterator it=it_begin; it<it_end; it++)
	{
		it_begin2 = lower_bound( qcal_tray.begin(), qcal_tray.end(), HKL_Q(NRVec<Int4>(), *it*(1.0 - resol) ) );
		it_end2 = upper_bound( it_begin2, qcal_tray.end(), HKL_Q(NRVec<Int4>(), *it*(1.0 + resol) ) );
		if( it_begin2 >= it_end2 ) return it;
		else closest_qcal_tray.push_back( closest_data(it_begin2, it_end2, *it) );
	}
	return it_end;
}

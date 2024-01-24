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
#ifndef _MATRIX_NBYN_HH_
#define _MATRIX_NBYN_HH_

#include"../RietveldAnalysisTypes.hh"
#include"../utility_data_structure/SymMat.hh"
#include"../utility_data_structure/nrutil_nr.hh"

inline void put_complement_set3(const Int4& i, const Int4& j, Int4& k)
{
	static const Int4 dim = 3;

	assert( i != j );
	assert( 0 <= i && i < dim );
	assert( 0 <= j && j < dim );

	bool flag_tray[dim] = {true, true, true};
	flag_tray[i]=false;
	flag_tray[j]=false;

	for(Int4 n=0; n<dim; n++)
	{
		if( flag_tray[n] )
		{
			k = n;
			return;
		}
	}
	k = -1; 	// To remove a warning.
}


inline void put_complement_set4(const Int4& i, const Int4& j, Int4& k, Int4& l)
{
	static const Int4 dim = 4;

	assert( i != j );
	assert( 0 <= i && i < dim );
	assert( 0 <= j && j < dim );

	bool flag_tray[dim] = {true, true, true, true};
	flag_tray[i]=false;
	flag_tray[j]=false;

	Int4 tray[2];
	for(Int4 n=0, index=0; index<2; n++)
	{
		if( flag_tray[n] ) tray[index++] = n;
	}
	k = tray[0];
	l = tray[1];
}


// 0<=i,j<ISIZE, i!=j.
inline NRMat<Int4> put_permutation_matrix(const Int4& ISIZE, const Int4& i, const Int4& j)
{
	NRMat<Int4> transmat(ISIZE,ISIZE,0);
	for(Int4 k=0; k<ISIZE; k++) transmat[k][k] = 1;
	transmat[i][i] = 0;
	transmat[j][j] = 0;
	transmat[i][j] = 1;
	transmat[j][i] = 1;

	return transmat;
}
// 0<=i<4, 0<=j<4, i!=j.
inline const NRMat<Int4>& put_permutation_matrix_dim_4(const Int4& i, const Int4& j)
{
	const NRMat<Int4> *trans_mat2;

// 2013.4.3 VIC Tamura. (For programs compiled with Visual Studio C++)
#ifdef _OPENMP
#pragma omp critical
#endif
    {
    	static const NRMat<Int4> trans_mat2_[6]
    	  = {
    		  put_permutation_matrix(4, 0, 1),
    		  put_permutation_matrix(4, 0, 2),
    		  put_permutation_matrix(4, 0, 3),
    		  put_permutation_matrix(4, 1, 2),
    		  put_permutation_matrix(4, 1, 3),
    		  put_permutation_matrix(4, 2, 3)
    	  };
          trans_mat2 = trans_mat2_;
    }

	assert( i < j );

	return trans_mat2[ i*(5-i)/2 + j-1 ];
}


// 0<=i,j<ISIZE, i!=j.
static NRMat<Int4> put_reduct_mat(const Int4& ISIZE, const Int4& i, const Int4& j)
{
	assert( 3 <= ISIZE && ISIZE <= 4 );
	if( ISIZE == 3 )
	{
		Int4 k;
		put_complement_set3(i, j, k);

		NRMat<Int4> transmat(3,3,0);
		transmat[i][i] = -1;
		transmat[j][j] = 1;
		transmat[k][k] = 1; transmat[k][i] = 2;

		return transmat;
	}
	else
	{
		Int4 k, l;
		put_complement_set4(i, j, k, l);

		NRMat<Int4> transmat(4,4,0);
		transmat[i][i] = -1;
		transmat[j][j] = 1;
		transmat[k][k] = 1; transmat[k][i] = 1;
		transmat[l][l] = 1; transmat[l][i] = 1;

		return transmat;
	}
}


inline const NRMat<Int4>& put_reduction_matrix(const Int4& ISIZE, const Int4& i, const Int4& j)
{
	assert( i < j );
	assert( 3 <= ISIZE && ISIZE <= 4 );

	// 2014.11.6 VIC Tamura. (For programs compiled with Visual Studio C++)
	const NRMat<Int4>* trans_mat3;
	const NRMat<Int4>* trans_mat4;
#ifdef _OPENMP
#pragma omp critical
#endif
{
    static const NRMat<Int4> trans_mat3_[3]
      = {
        put_reduct_mat(3, 0, 1),
        put_reduct_mat(3, 0, 2),
        put_reduct_mat(3, 1, 2),
    };
    static const NRMat<Int4> trans_mat4_[6]
      = {
        put_reduct_mat(4, 0, 1),
        put_reduct_mat(4, 0, 2),
        put_reduct_mat(4, 0, 3),
        put_reduct_mat(4, 1, 2),
        put_reduct_mat(4, 1, 3),
        put_reduct_mat(4, 2, 3)
    };
    trans_mat3 = trans_mat3_;
    trans_mat4 = trans_mat4_;
}

	if( ISIZE == 3 ) return trans_mat3[ i*(3-i)/2 + j-1 ];
	else return trans_mat4[ i*(5-i)/2 + j-1 ];
}


inline NRMat<Int4> put_transform_matrix_rowNtoNplus1(const Int4& isize)
{
	assert(isize >= 0);
	NRMat<Int4> TransMat(isize+1,isize,0);
	for(Int4 k=0; k<isize; k++) TransMat[k][k] = 1;
	for(Int4 k=0; k<isize; k++) TransMat[isize][k] = -1;

	return TransMat;
}


inline const NRMat<Int4>& put_transform_matrix_row3to4()
{
	static const NRMat<Int4> mat43 = put_transform_matrix_rowNtoNplus1(3);
	return mat43;
}


inline NRMat<Int4> put_transform_matrix_row3to4(const NRMat<Int4>& rhs)
{
	assert( rhs.nrows() == 3 );
	
	const Int4 icol = rhs.ncols();
	
	NRMat<Int4> TransMat(4, icol);
	for(Int4 k=0; k<3; k++)
	{
		for(Int4 k2=0; k2<icol; k2++)
		{
			TransMat[k][k2] = rhs[k][k2];
		}
	}
	for(Int4 k2=0; k2<icol; k2++)
	{
		TransMat[3][k2] = -( rhs[0][k2] + rhs[1][k2] + rhs[2][k2] );
	}

	return TransMat;
}


inline NRMat<Int4> put_transform_matrix_row4to3(const NRMat<Int4>& rhs)
{
	assert( rhs.nrows() == 4 );
	
	const Int4 icol = rhs.ncols();
	NRMat<Int4> TransMat(3, icol);
	for(Int4 k=0; k<3; k++)
	{
		for(Int4 k2=0; k2<icol; k2++)
		{
			TransMat[k][k2] = rhs[k][k2];
		}
	}
	return TransMat;
}


template<class T>
inline SymMat<T> put_sym_matrix_sizeNplus1toN(const SymMat<T>& rhs)
{
	assert( rhs.size() > 0 );

	static const Int4 isize = rhs.size() - 1;
	SymMat<T> ans(isize);
	
	for(Int4 k=0; k<isize; k++)
	{
		for(Int4 k2=k; k2<isize; k2++)
		{
			ans(k,k2) = rhs(k,k2);
		}
	}
	return ans;
}

template<class T>
inline SymMat<T> put_sym_matrix_sizeNtoNplus1(const SymMat<T>& rhs)
{
	static const Int4 isize = rhs.size();
	SymMat<T> ans(isize+1);
	
	// Copy.
	for(Int4 k=0; k<isize; k++)
	{
		for(Int4 k2=k; k2<isize; k2++)
		{
			ans(k,k2) = rhs(k,k2);
		}
	}

	ans(isize,isize) = 0.0;
	for(Int4 k=0; k<isize; k++)
	{
		ans(k, isize) = 0.0;
		for(Int4 k2=0; k2<isize; k2++)
		{
			ans(k, isize) -= rhs(k, k2);
		}
		ans(isize,isize) -= ans(k, isize);
	}

	return ans;
}

#endif

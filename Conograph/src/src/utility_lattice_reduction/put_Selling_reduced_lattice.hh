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
#ifndef _PUT_SELLING_REDUCED_LATTICE_HH_
#define _PUT_SELLING_REDUCED_LATTICE_HH_

#include "matrix_NbyN.hh"
#include "../utility_data_structure/index_set.hh"
#include "../utility_func/transform_sym_matrix.hh"

// Transform a symmetric matrix S such that the diagonals are sorted into ascending order.
template<class T>
inline void putMatrixToMoveSmallerDiagonalLeftUpper(const SymMat<T>& S, NRMat<Int4>& trans_mat2)
{
	const Int4 ISIZE = S.size();
	vector< index_set<T> > tray(ISIZE);
	for(Int4 k=0; k<ISIZE; k++)
	{
		tray[k].index = k;
		tray[k].element = S(k, k);
	}
	sort( tray.begin(), tray.end() );
	
	trans_mat2 = NRMat<Int4>(ISIZE, ISIZE, 0);
	for(Int4 k=0; k<ISIZE; k++) trans_mat2[k][ tray[k].index ] = 1;
}


template<class T>
inline void putMatrixToMoveLargerDiagonalLeftUpper(const SymMat<T>& S, NRMat<Int4>& trans_mat2)
{
	const Int4 ISIZE = S.size();
	vector< index_set<T> > tray(ISIZE);
	for(Int4 k=0; k<ISIZE; k++)
	{
		tray[k].index = k;
		tray[k].element = S(k, k);
	}
	sort( tray.begin(), tray.end() );
	
	trans_mat2 = NRMat<Int4>(ISIZE, ISIZE, 0);
	for(Int4 k=0, k2=ISIZE-1; k<ISIZE; k++, k2--) trans_mat2[k2][ tray[k].index ] = 1;
}


template<class T>
inline void moveSmallerDiagonalLeftUpper(SymMat<T>& S)
{
	NRMat<Int4> trans_mat2;
	putMatrixToMoveSmallerDiagonalLeftUpper(S, trans_mat2);
	S = transform_sym_matrix(trans_mat2, S);
}


template<class T>
inline void moveSmallerDiagonalLeftUpper(SymMat<T>& S, NRMat<Int4>& trans_mat)
{
	NRMat<Int4> trans_mat2;
	putMatrixToMoveSmallerDiagonalLeftUpper(S, trans_mat2);
	S = transform_sym_matrix(trans_mat2, S);
	trans_mat =  mprod(trans_mat2, trans_mat);
}


template<class T>
inline void moveLargerDiagonalLeftUpper(SymMat<T>& S, NRMat<Int4>& trans_mat)
{
	NRMat<Int4> trans_mat2;
	putMatrixToMoveLargerDiagonalLeftUpper(S, trans_mat2);
	S = transform_sym_matrix(trans_mat2, S);
	trans_mat =  mprod(trans_mat2, trans_mat);
}


// Transform S by V such that the nondiagonal elements of V * S * transpose(V) are all not positive.
// On output, trans_mat and inv_trans_mat vector_ortho_trans_mat are replaced
// by V * transmat, inv_trans_mat = inv_tran_mat * V^-1, vector_ortho_trans_mat * V^-1.
// S is the matrix whose nondiagonal elements are all not positive.
template<class T>
bool put_Selling_reduced_dim_less_than_4(SymMat<T>& S, NRMat<Int4>& trans_mat, const T Min_Sii = 0)
{
	static const Int4 MAXIT = 10000;
	const Int4 ISIZE = S.size();
	assert( ISIZE <= 4 );
	static const T zerro = 0;
	
	Int4 itnum = 0;
	Int4 i, j;
	NRMat<Int4> trans_mat2(ISIZE,ISIZE);
	while(itnum < MAXIT)
	{
		itnum++;
		
		for(i=0; i<ISIZE; i++)
			if( S(i,i) < Min_Sii ) return false;
		
		for(i=0; i<ISIZE; i++)
		{
			for(j=i+1; j<ISIZE; j++)
			{
				if( zerro < S(i,j) ) break;
			}
			if( j < ISIZE ) break;
		}
		
		if( i >= ISIZE )
		{
			return true; // S is Selling-reduced.
		}

		trans_mat2 = put_reduction_matrix(S.size(), i, j);
		S = transform_sym_matrix(trans_mat2, S);

		trans_mat = mprod(trans_mat2, trans_mat);
	}
	
	return false;
}


template<class T>
inline bool put_Selling_reduced_dim_less_than_4(const SymMat<T>& S,
										SymMat<T>& S_super, NRMat<Int4>& TransMat, const T Min_Sii = 0)
{
	TransMat = put_transform_matrix_rowNtoNplus1(S.size());
	
	// Copy.
	S_super = transform_sym_matrix(TransMat, S);

	return put_Selling_reduced_dim_less_than_4(S_super, TransMat, Min_Sii);
}


template<class T>
inline bool put_Selling_reduced_matrix(const SymMat<T>& S, SymMat<T>& S_super)
{
	NRMat<Int4> TransMat(S.size()+1, S.size());
	return put_Selling_reduced_dim_less_than_4(S, S_super, TransMat);
}

#endif /*SUPER_BASIS3_HH_*/

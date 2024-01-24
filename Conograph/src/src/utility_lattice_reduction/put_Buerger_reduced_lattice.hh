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
#ifndef PUT_Buerger_REDUCED_LATTICE_HH_
#define PUT_Buerger_REDUCED_LATTICE_HH_

#include "../utility_data_structure/SymMat.hh"
#include "../utility_data_structure/VCData.hh"
#include "../utility_data_structure/nrutil_nr.hh"
#include "../utility_lattice_reduction/put_Selling_reduced_lattice.hh"
#include "../point_group/enumPointGroup.hh"
#include "../point_group/point_gp_data.hh"
#include "../centring_type/enumCentringType.hh"
#include "../bravais_type/BravaisType.hh"


template<class T>
static void arrangeNondiagonalSign(SymMat<T>& inv_S_red, NRMat<Int4>& trans_mat)
{
	static const T zerro = 0;

	NRMat<Int4> mat2(3,3,0);
	mat2[0][0] = 1;
	mat2[1][1] = 1;
	mat2[2][2] = 1;
	
	if( zerro < inv_S_red(0,1) )
	{
		if( zerro < inv_S_red(0,2) )
		{
			if( zerro < inv_S_red(1,2) ) return;
			else // if( inv_S_red(1,2) <= zerro )
			{
				mat2[0][0] = -1;
				inv_S_red(0,1) *= -1;
				inv_S_red(0,2) *= -1;
			}
		}
		else // if( inv_S_red(0,2) <= zerro )
		{
			if( zerro <= inv_S_red(1,2) )
			{
				mat2[1][1] = -1;
				inv_S_red(0,1) *= -1;
				inv_S_red(1,2) *= -1;
			}
			else // if( inv_S_red(1,2) < zerro )
			{
				if( zerro <= inv_S_red(0,2) )  // zerro == inv_S_red(0,2)
				{
					mat2[0][0] = -1;
					inv_S_red(0,1) *= -1;
				}
				else // inv_S_red(0,2) < zerro
				{
					mat2[2][2] = -1;
					inv_S_red(0,2) *= -1;
					inv_S_red(1,2) *= -1;
				}
			}
		}
	}
	else //if( inv_S_red(0,1) <= zerro )
	{
		if( zerro <= inv_S_red(0,2) )
		{
			if( zerro <= inv_S_red(1,2) )
			{
				mat2[2][2] = -1;
				inv_S_red(0,2) *= -1;
				inv_S_red(1,2) *= -1;
			}
			else // if( inv_S_red(1,2) < zerro )
			{
				if( inv_S_red(0,2) <= zerro ) return;
				else if( zerro <= inv_S_red(0,1) )	// zerro == inv_S_red(0,1).
				{
					mat2[0][0] = -1;
					inv_S_red(0,2) *= -1;
				}
				else // inv_S_red(0,1) < zerro && inv_S_red(0,2) > zerro.
				{
					mat2[1][1] = -1;
					inv_S_red(0,1) *= -1;
					inv_S_red(1,2) *= -1;
				}
			}
		}
		else // if( inv_S_red(0,2) < zerro )
		{
			if( zerro < inv_S_red(1,2) )
			{
				if( zerro <= inv_S_red(0,1) )	// inv_S_red(0,1) = zerro.
				{
					mat2[1][1] = -1;
					inv_S_red(1,2) *= -1;
				}
				else	// inv_S_red(0,1) < zerro.
				{
					mat2[0][0] = -1;
					inv_S_red(0,1) *= -1;
					inv_S_red(0,2) *= -1;
				}
			}
			else // if( inv_S_red(1,2) <= zerro )
				return;
		}
	}
	trans_mat = mprod(mat2, trans_mat);
}


// inv_S_red = trans_mat2*put_transform_matrix_34() * inv_S_super * Transpose(trans_mat2*put_transform_matrix_34()).
template<class T>
void putBuergerReducedMatrix(const SymMat<T>& inv_S_super, const bool& inv_flag, SymMat<T>& inv_S_red,
		NRMat<Int4>& trans_mat2)
{
	static const T zerro = 0;

	assert( inv_S_super.size() == 4 );
	assert( inv_S_red.size() == 3 );
	assert( inv_S_super(0,0) <= inv_S_super(1,1) );
	assert( inv_S_super(1,1) <= inv_S_super(2,2) );
	assert( inv_S_super(2,2) <= inv_S_super(3,3) );
	assert( inv_S_super(0,1) <= zerro
			&& inv_S_super(0,2) <= zerro
			&& inv_S_super(0,3) <= zerro
			&& inv_S_super(1,2) <= zerro
			&& inv_S_super(1,3) <= zerro
			&& inv_S_super(2,3) <= zerro );
	
	trans_mat2 = NRMat<Int4>(3,3,0);
	trans_mat2[0][0] = 1;
	trans_mat2[1][1] = 1;
	trans_mat2[2][2] = 1;
	
	if( inv_S_super(0,0) + inv_S_super(0,1)*2 < zerro )
	{
		trans_mat2[1][0] = 1;
	}
	else if( inv_S_super(0,0) + inv_S_super(0,2)*2 < zerro )
	{
		trans_mat2[2][0] = 1;
	}
	else if( inv_S_super(1,1) + inv_S_super(1,2)*2 < zerro )
	{
		trans_mat2[2][1] = 1;
	}
	else if( inv_S_super(0,0) + inv_S_super(1,1) + (inv_S_super(0,1) + inv_S_super(0,2) + inv_S_super(1,2) )*2 < zerro )
	{
		trans_mat2[2][0] = -1;
		trans_mat2[2][1] = -1;
		trans_mat2[2][2] = -1;
	}
	inv_S_red = transform_sym_matrix(trans_mat2, put_sym_matrix_sizeNplus1toN(inv_S_super));
	if( inv_flag ) moveSmallerDiagonalLeftUpper(inv_S_red, trans_mat2);
	else moveLargerDiagonalLeftUpper(inv_S_red, trans_mat2);
	arrangeNondiagonalSign(inv_S_red, trans_mat2);
}



template <class T>
void cal_average_crystal_system(const ePointGroup& epg, SymMat<T>& ans)
{
	if(epg == Ci) return;
	else if(epg == C2h_X)
	{
		ans(0,1) = 0;
		ans(0,2) = 0;
	}
	else if(epg == C2h_Y)
	{
		ans(0,1) = 0;
		ans(1,2) = 0;
	}
	else if(epg == C2h_Z)
	{
		ans(0,2) = 0;
		ans(1,2) = 0;
	}
	else if(epg == D2h)
	{
		ans(0,1) = 0;
		ans(0,2) = 0;
		ans(1,2) = 0;
	}
	else if(epg == D4h_X)
	{
		ans(0,1) = 0;
		ans(0,2) = 0;
		ans(1,2) = 0;
		ans(1,1) = (ans(1,1)+ans(2,2))/2;
		ans(2,2) = ans(1,1);
	}
	else if(epg == D4h_Y)
	{
		ans(0,1) = 0;
		ans(0,2) = 0;
		ans(1,2) = 0;
		ans(0,0) = (ans(0,0)+ans(2,2))/2;
		ans(2,2) = ans(0,0);
	}
	else if(epg == D4h_Z)
	{
		ans(0,1) = 0;
		ans(0,2) = 0;
		ans(1,2) = 0;
		ans(0,0) = (ans(0,0)+ans(1,1))/2;
		ans(1,1) = ans(0,0);
	}
	else if(epg == D31d_rho)
	{
		ans(0,0) = (ans(0,0)+ans(1,1)+ans(2,2))/3;
		ans(1,1) = ans(0,0);
		ans(2,2) = ans(0,0);
		ans(0,1) = (ans(0,1)+ans(0,2)+ans(1,2))/3;
		ans(0,2) = ans(0,1);
		ans(1,2) = ans(0,1);
	}
	else if(epg == D3d_1_hex || epg == D6h)
	{
		ans(0,0) = (ans(0,0)+ans(1,1))/2;
		ans(1,1) = ans(0,0);
		ans(0,1) = ans(0,0)/2;
		ans(0,2) = 0;
		ans(1,2) = 0;
	}
	else if(epg == Oh)
	{
		ans(0,0) = (ans(0,0)+ans(1,1)+ans(2,2))/3;
		ans(1,1) = ans(0,0);
		ans(2,2) = ans(0,0);
		ans(0,1) = 0;
		ans(0,2) = 0;
		ans(1,2) = 0;
	}
	else
	{
		assert( false );
	}
};

inline void putBuergerReducedMatrix(const SymMat<Double>& inv_S_super, SymMat<Double>& inv_S_red, NRMat<Int4>& trans_mat2)
{
	putBuergerReducedMatrix(inv_S_super, true, inv_S_red, trans_mat2);

#ifdef DEBUG
	assert( ( inv_S_red(0,1) <= 0.0 && inv_S_red(0,2) <= 0.0 && inv_S_red(1,2) <= 0.0 )
			|| ( 0.0 < inv_S_red(0,1) && 0.0 < inv_S_red(0,2) && 0.0 < inv_S_red(1,2) ) );
	assert( inv_S_red(1,1)*0.9999 < inv_S_red(2,2) );
	assert( inv_S_red(0,0)*0.9999 < inv_S_red(1,1) );
	assert( inv_S_red(0,1) * (-1.9999) < inv_S_red(0,0)
					&& inv_S_red(0,1) * 1.9999 < inv_S_red(0,0)
					&& inv_S_red(0,2) * (-1.9999) < inv_S_red(0,0)
					&& inv_S_red(0,2) * 1.9999 < inv_S_red(0,0)
					&& inv_S_red(1,2) * (-1.9999) < inv_S_red(1,1)
					&& inv_S_red(1,2) * 1.9999 < inv_S_red(1,1) );
	assert( 0.0 < inv_S_red(0,0) + inv_S_red(1,1) + ( inv_S_red(0,1) + inv_S_red(0,2) + inv_S_red(1,2) ) * 1.9999
					&& 0.0 < inv_S_red(0,0) + inv_S_red(1,1) + ( inv_S_red(0,1) - inv_S_red(0,2) - inv_S_red(1,2) ) * 1.9999
					&& 0.0 < inv_S_red(0,0) + inv_S_red(1,1) - ( inv_S_red(0,1) - inv_S_red(0,2) + inv_S_red(1,2) ) * 1.9999
					&& 0.0 < inv_S_red(0,0) + inv_S_red(1,1) - ( inv_S_red(0,1) + inv_S_red(0,2) - inv_S_red(1,2) ) * 1.9999 );
#endif
}


inline void putBuergerReducedMatrix(const SymMat<VCData>& S_super, SymMat<VCData>& S_red, NRMat<Int4>& trans_mat2)
{
	putBuergerReducedMatrix(S_super, false, S_red, trans_mat2);

#ifdef DEBUG
	assert( ( S_red(0,1).Value() <= 0.0 && S_red(0,2).Value() <= 0.0 && S_red(1,2).Value() <= 0.0 )
			|| ( 0.0 < S_red(0,1).Value() && 0.0 < S_red(0,2).Value() && 0.0 < S_red(1,2).Value() ) );
	assert( S_red(1,1).Value()*0.9999 < S_red(0,0).Value() );
	assert( S_red(2,2).Value()*0.9999 < S_red(1,1).Value() );
	assert( S_red(0,1).Value() * (-1.9999) < S_red(1,1).Value()
					&& S_red(0,1).Value() * 1.9999 < S_red(1,1).Value()
					&& S_red(0,2).Value() * (-1.9999) < S_red(2,2).Value()
					&& S_red(0,2).Value() * 1.9999 < S_red(2,2).Value()
					&& S_red(1,2).Value() * (-1.9999) < S_red(2,2).Value()
					&& S_red(1,2).Value() * 1.9999 < S_red(2,2).Value() );
	assert( 0.0 < S_red(1,1).Value() + S_red(2,2).Value() + ( S_red(0,1) + S_red(0,2) + S_red(1,2) ).Value() * 1.9999
					&& 0.0 < S_red(1,1).Value() + S_red(2,2).Value() + ( S_red(0,1) - S_red(0,2) - S_red(1,2) ).Value() * 1.9999
					&& 0.0 < S_red(1,1).Value() + S_red(2,2).Value() - ( S_red(0,1) - S_red(0,2) + S_red(1,2) ).Value() * 1.9999
					&& 0.0 < S_red(1,1).Value() + S_red(2,2).Value() - ( S_red(0,1) + S_red(0,2) - S_red(1,2) ).Value() * 1.9999 );
#endif
}


template<class T>
void putBuergerReducedMonoclinicP(const Int4& i, const Int4& j,
		SymMat<T>& S_red, NRMat<Int4>& trans_mat2)
{
	static const T zerro = 0;

	assert(S_red.size()==3);
	assert(trans_mat2.ncols()==3);
	assert(i < j);
	const Int4 irow = trans_mat2.nrows();

	if( S_red(i,j) < zerro )
	{
		S_red(i,j) *= -1;
		for(Int4 l=0; l<irow; l++)
		{
			trans_mat2[l][i] *= -1;
		}
	}

	do{
		if( S_red(j,j) < S_red(i,j) * 2 )
		{
			// S_red(i,j).Value() <= S_red(j,j).Value() * m
			// i :  -1  m  0
			// j :   0  1  0
			//   :   0  0  1
			const Int8 m = iceil( S_red(i,j) / S_red(j,j) );

			S_red(i,i) += ( S_red(j,j) * m - S_red(i,j) * 2 ) * m;
			assert( zerro < S_red(i,i) );
			S_red(i,j) = S_red(j,j) * m - S_red(i,j);

			for(Int4 l=0; l<irow; l++)
			{
				trans_mat2[l][j] += trans_mat2[l][i]*m;
			}
			for(Int4 l=0; l<irow; l++)
			{
				trans_mat2[l][i] *= -1;
			}
		}

		if( S_red(i,i) < S_red(i,j) * 2 )
		{
			// S_red(i,j).Value() <= S_red(i,i).Value() * n
			// i : 1  0  0
			// j : n -1  0
			//   : 0  0  1
			const Int8 n = iceil( S_red(i,j) / S_red(i,i) );

			S_red(j,j) += ( S_red(i,i) * n - S_red(i,j) * 2 ) * n;
			assert( zerro < S_red(j,j) );
			S_red(i,j) = S_red(i,i) * n - S_red(i,j);

			for(Int4 l=0; l<irow; l++)
			{
				trans_mat2[l][i] += trans_mat2[l][j]*n;
			}
			for(Int4 l=0; l<irow; l++)
			{
				trans_mat2[l][j] *= -1;
			}
		}
	}
	while( S_red(j,j) < S_red(i,j) * 2 || S_red(i,i) < S_red(i,j) * 2 );

	if( S_red(i,i) < S_red(j,j) )
	{
		const Int4 k = put_complement_set3(i, j);
		swap(S_red(i,i), S_red(j,j));
		swap(S_red(i,k), S_red(j,k));
		for(Int4 l=0; l<irow; l++)
		{
			swap(trans_mat2[l][i], trans_mat2[l][j]);
		}
	}
}


inline Double operator/(const VCData& lhs, const VCData& rhs)
{
	assert(rhs.Value() != 0.0);
	return lhs.Value() / rhs.Value();
}

template<class T>
void putBuergerReducedMonoclinicB(
		const BravaisType& monoclinic_b_type,
		SymMat<T>& S_red, NRMat<Int4>& trans_mat2)
{
	const Int4 ibase_axis = monoclinic_b_type.enumBASEaxis();
	const Int4 iabc_axis = monoclinic_b_type.enumABCaxis();
	const Int4 ir = put_complement_set3(iabc_axis, ibase_axis);

	static const T zerro = 0;

	assert(S_red.size()==3);
	assert(trans_mat2.nrows()==0 || trans_mat2.ncols()==3);
	const Int4 irow = trans_mat2.nrows();

	if( S_red( ibase_axis, ir ) < zerro )
	{
		S_red( ibase_axis, ir ) *= -1;
		for(Int4 l=0; l<irow; l++)
		{
			trans_mat2[l][ir] *= -1;
		}
	}

	do{
		if( S_red( ir, ir ) < S_red( ibase_axis, ir ) )
		{
			// S_red( ibase_axis, ir ).Value() <= S_red( ir, ir ).Value() * m2
			//         ir :   1  0  0
			// ibase_axis :  m2 -1  0
			//            :   0  0  1
			const Int8 m1 = iceil( ( S_red( ibase_axis, ir ) / S_red( ir, ir ) ) * 0.5 );
			const Int8 m2 = m1*2;

			S_red( ibase_axis, ibase_axis ) += ( S_red( ir, ir ) * m2 - S_red( ibase_axis, ir ) * 2 ) * m2;
			assert( zerro < S_red( ibase_axis, ibase_axis ) );
			S_red( ibase_axis, ir ) = S_red( ir, ir ) * m2 - S_red( ibase_axis, ir );

			for(Int4 l=0; l<irow; l++)
			{
				trans_mat2[l][ir] += trans_mat2[l][ibase_axis]*m2;
			}
			for(Int4 l=0; l<irow; l++)
			{
				trans_mat2[l][ibase_axis] *= -1;
			}
		}

		if( S_red( ibase_axis, ibase_axis ) < S_red( ibase_axis, ir ) * 2  )
		{
			// S_red( ibase_axis, ir ).Value() <= S_red( ibase_axis, ibase_axis ).Value() * n
			//         ir :-1  n  0
			// ibase_axis : 0  1  0
			//            : 0  0  1
			const Int4 n = iceil( S_red( ibase_axis, ir ) / S_red( ibase_axis, ibase_axis ) );

			S_red( ir, ir ) += ( S_red( ibase_axis, ibase_axis ) * n - S_red( ibase_axis, ir ) * 2 ) * n;
			assert( zerro < S_red( ir, ir ) );
			S_red( ibase_axis, ir ) = S_red( ibase_axis, ibase_axis ) * n - S_red( ibase_axis, ir );

			for(Int4 l=0; l<irow; l++)
			{
				trans_mat2[l][ibase_axis] += trans_mat2[l][ir]*n;
			}
			for(Int4 l=0; l<irow; l++)
			{
				trans_mat2[l][ir] *= -1;
			}
		}
	}
	while( S_red( ir, ir ) < S_red( ibase_axis, ir ) || S_red( ibase_axis, ibase_axis ) < S_red( ibase_axis, ir ) * 2 );
}


template<class T>
void putBuergerReducedOrthorhombic(const eCentringType& brat,
		SymMat<T>& S_red, NRMat<Int4>& trans_mat)
{
	assert( brat != BaseX );
	assert( brat != BaseY );

	if( brat == BaseZ )
	{
		if( S_red(0,0) < S_red(1,1) )
		{
			S_red = transform_sym_matrix(put_matrix_YXZ(), S_red);
			trans_mat = mprod( trans_mat, put_matrix_YXZ() );
		}
	}
	else
	{
		NRMat<Int4> trans_mat2 = identity_matrix<Int4>(3);
		moveLargerDiagonalLeftUpper(S_red, trans_mat2);
		trans_mat = mprod(trans_mat, transpose(trans_mat2));	// inverse(trans_mat2) = transpose(trans_mat2).
	}
}

#endif /*PUT_MINKOWSKI_REDUCED_LATTICE_HH_*/

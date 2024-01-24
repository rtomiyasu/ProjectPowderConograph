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
#ifndef SymmetricOperation_HH_
#define SymmetricOperation_HH_

#include "../RietveldAnalysisTypes.hh"
#include "../point_group/enumPointGroup.hh"
#include "../utility_data_structure/nrutil_nr.hh"
#include "enumSymmetricOperation.hh"

static const Int4 ISIZE = 4;

typedef struct{
	UChar Perm[ISIZE];
	Char Sign[ISIZE];
} SymmetricOperation;


inline void putMatrixForm(const SymmetricOperation& sop, NRMat<Int4>& mat)
{
	mat = NRMat<Int4>(3,3,0);
	for(Int4 i=0; i<3; i++)
		if(sop.Perm[i]!=3) mat[sop.Perm[i]][i] = sop.Sign[i];
	if(sop.Perm[3]!=3)
	{
		mat[sop.Perm[3]][0] = -sop.Sign[3];
		mat[sop.Perm[3]][1] = -sop.Sign[3];
	}
}

inline SymmetricOperation operator*(const SymmetricOperation& lhs, const SymmetricOperation& rhs)
{
	SymmetricOperation ans;
	for(Int4 i=0; i<ISIZE; i++){
		ans.Perm[i]=lhs.Perm[ rhs.Perm[i] ];
		ans.Sign[i]=lhs.Sign[ rhs.Perm[i] ] * rhs.Sign[i];
	}
	return ans;
}

inline SymmetricOperation Inverse(const SymmetricOperation& rhs)
{
	SymmetricOperation ans;
	for(Int4 i=0; i<ISIZE; i++){
		ans.Perm[ rhs.Perm[i] ]=i;
		ans.Sign[ rhs.Perm[i] ]=rhs.Sign[i];
	} 
	return ans;
}


inline bool operator==(const SymmetricOperation& lhs, const SymmetricOperation& rhs)
{
	for(Int4 i=0; i<3; i++)
	{
		if( lhs.Perm[i]!=rhs.Perm[i] ) return false;
		if( lhs.Sign[i]!=rhs.Sign[i] ) return false;
	} 
	if(lhs.Perm[3]!=3) 
	{
		if( lhs.Perm[3]!=rhs.Perm[3] ) return false;
		if( lhs.Sign[3]!=rhs.Sign[3] ) return false;
	}
	return true;
}


const SymmetricOperation& change_enum_to_data(const eSymmetricOperation&);
bool change_data_to_enum(const SymmetricOperation&, eSymmetricOperation&);

const string& Name(const eSymmetricOperation& num);
const ePointGroup& generateGroup(const eSymmetricOperation&);
inline const ePointGroup& generateGroup(const SymmetricOperation& sop)
{
	eSymmetricOperation esym;
	if( !change_data_to_enum(sop, esym) )
		if( !change_data_to_enum(sop*sop, esym) ) assert( false );
	return generateGroup(esym);
}

void putAllElement(const ePointGroup& epg, vector<SymmetricOperation>&);

void putSumAllElement(const ePointGroup& epg, NRMat<Int4>& ans);

const SymmetricOperation& E();
const SymmetricOperation& I();

#endif /*SymmetricOperation_HH_*/

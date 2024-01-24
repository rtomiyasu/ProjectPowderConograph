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
#include <algorithm>
#include "PGNormalSeriesTray.hh"
#include"../point_group/coset_representative_data.hh"

// On output, the size of coset_rep is ilength - 1.
static void CosetRepresentativeOfNormalSeries(const ePointGroup& epg, Int4& ilength,
vector<Int4>& index, vector<eSymmetricOperation>& coset_rep)
{
	index.clear();
	coset_rep.clear();
	
	// Set normal series from epg to C1.
	ePointGroup epg2 = epg;
	eGroupToMaxSubgp egp_subgp;
	ilength = 1;
	while(epg2 != C1)
	{
		egp_subgp = enumMaxNormalSubgroup(epg2);

		index.resize(ilength);
		coset_rep.resize(ilength);
		
		CosetRepresentativeMaxSubgp(egp_subgp, index[ilength-1], coset_rep[ilength-1]);
		ilength++;

		epg2 = enumLowerGroup( egp_subgp );
	}
}


PGNormalSeriesTray::PGNormalSeriesTray(const ePointGroup& num) : m_epg(num)
{
	CosetRepresentativeOfNormalSeries(m_epg, m_length, m_index, m_coset_rep);  // Normal series of m_sg.
}

PGNormalSeriesTray::~PGNormalSeriesTray()
{
}


Int4 PGNormalSeriesTray::HKLEquiv(const MillerIndex3& hkl,
const Int4& ilength, 
const vector<Int4>::const_iterator& it_index,
const vector<SymmetricOperation>::const_iterator& it_coset_rep,
vector<MillerIndex3>& hkl_equiv) const
{
	if(ilength>1)
	{
		const Int4 multifactor = HKLEquiv(hkl, ilength-1, it_index+1, it_coset_rep+1, hkl_equiv);
		const Int4 isize = hkl_equiv.size();
		
		MillerIndex3 hkl2 = (*it_coset_rep)*hkl;
		Int4 j = distance( hkl_equiv.begin(), find( hkl_equiv.begin(), hkl_equiv.end(), hkl2 ) );

		if( j<isize )
		{
			return multifactor*(*it_index);
		}
		for(Int4 i=0; i<*it_index-1; i++)
		{
			vector<MillerIndex3> hkl_equiv2;
			HKLEquiv(hkl2, ilength-1, it_index+1, it_coset_rep+1, hkl_equiv2);
			hkl_equiv.insert( hkl_equiv.end(), hkl_equiv2.begin(), hkl_equiv2.end() );
			hkl2 = (*it_coset_rep)*hkl2;
		}
		return multifactor;
	}
	else
	{
		hkl_equiv.clear();
		hkl_equiv.resize(1, hkl);
		return 1; 
	}
}

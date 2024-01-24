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
#ifdef _OPENMP
 # include <omp.h>
#endif
#include <assert.h>
#include "zerror_type/error_out.hh"
#include "GenerateRelation.hh"
#include "utility_data_structure/Bud.hh"
#include "utility_data_structure/VCData.hh"
#include "utility_func/stopx.hh"

GenerateRelation::GenerateRelation(const Int4& max_data_num) : m_max_data_num (max_data_num)
{	
	itype = -1;
}

GenerateRelation::~GenerateRelation()
{
}

void GenerateRelation::initialize_type2() 
{
	itype = 2;

	m_vecdat.clear();
	const vector<QData>& qdata = VCData::putPeakQData();

	const Int4 peak_num = min(m_max_data_num, (Int4)qdata.size());
	
	// Set Qi+Qj, 2*(Qi+Qj).
	for(Int4 k=0; k<peak_num; k++)
	{
		VCData vctray(k, 1);

		for(Int4 j=0; j<=k; j++)
		{
			VCData vctray2(j, 1);
			vctray2 += vctray;
			m_vecdat.push_back( vctray2 );
		}
	}

	stable_sort( m_vecdat.begin(), m_vecdat.end() );
}


void GenerateRelation::initialize_type3() 
{
	itype = 3;

	m_vecdat.clear();

	const Int4 peak_num = min(m_max_data_num, (Int4)VCData::putPeakQData().size());
	
	// Set Qi+3*Qj, 3*Qi+Qj.
	for(Int4 k=0; k<peak_num; k++)
	{
		VCData vctray(k, 1);

		for(Int4 j=0; j<k; j++)
		{
			VCData vctray2(j, 1);
			m_vecdat.push_back( vctray*3 + vctray2 );
			m_vecdat.push_back( vctray + vctray2*3 );
		}

		m_vecdat.push_back( vctray * 4 );
	}

	stable_sort( m_vecdat.begin(), m_vecdat.end() );
}


// On output, lhs = rhs, lhs = 2 *(Qi+Qj) and rhs = Qk+Ql.
void GenerateRelation::putRelationCandidate_type2(const Double& cv2,
		vector<Bud>& ans) const
{
	ans.clear();
	if( m_vecdat.empty() ) return;

	const Int4 max_size = distance(m_vecdat.begin(),
			upper_bound(m_vecdat.begin(), m_vecdat.end(), VCData::putVCData( (m_vecdat.rbegin()->Value()+sqrt(m_vecdat.rbegin()->Variance()*cv2))*0.5) ));

	Bud budex;

	assert( itype == 2 );

	try{

/*   2011.10.19 VIC Tamura */
Int4 LOOP_COUNTER = 0;

#ifdef _OPENMP
		#pragma omp for private(budex)
#endif
		for(Int4 i_=0; i_<max_size; i_++)
		{
/*   2011.10.19 VIC Tamura */
SET_PROGRESS(100*(LOOP_COUNTER++)/max_size, 2, 1); // critical, but works
if(IS_CANSELED()) continue;
			
			const vector<VCData>::const_iterator it=m_vecdat.begin()+i_;
			const Double it_error = sqrt(it->Variance()*cv2)*2.0;
			const vector<VCData>::const_iterator it2_begin=lower_bound(m_vecdat.begin(), m_vecdat.end(), VCData::putVCData(it->Value()*2.0-it_error));
			const vector<VCData>::const_iterator it2_end=upper_bound(it2_begin, m_vecdat.end(), VCData::putVCData(it->Value()*2.0+it_error));

			for(vector<VCData>::const_iterator it2 = it2_begin; it2<it2_end; it2++)
			{
				const VCData diff = (*it)*2 - *it2;
				if( diff.Value()*diff.Value() > cv2 * it2->Variance() ) continue;

				// lhs = 2 *(Qi+Qj)
				// rhs = Qk+Ql
				budex.setIndex( it->putVecCoef().begin()->first,
						        it->putVecCoef().rbegin()->first,
								it2->putVecCoef().begin()->first,
								it2->putVecCoef().rbegin()->first );
				if( budex.cross_product_312() <= 0.0
						|| budex.cross_product_412() <= 0.0 ) continue;

#ifdef _OPENMP
				#pragma omp critical
#endif
				{
					ans.push_back(budex);
				}
			}
		}
	}
	catch(bad_alloc& ball)
	{
		throw nerror(ball, __FILE__, __LINE__, __FUNCTION__);
	}

/*   2011.10.19 VIC Tamura */
CHECK_INTERRUPTION();
}



void GenerateRelation::putRelationCandidate_type3(const Double& cv2,
		vector<Bud>& ans) const
{
	const Int4 max_size = m_vecdat.size();

	ans.clear();
	Bud budex, budex2;
	set<Bud> found_bud_tray;

	static const type_coef i3 = VCData::putDenom()*3;
	
	assert( itype == 3 );

	try{

/*   2011.10.19 VIC Tamura */
Int4 LOOP_COUNTER = 0;

#ifdef _OPENMP
		#pragma omp parallel for private(budex, budex2)
#endif
	for(Int4 i_=0; i_<max_size; i_++)
	{
/*   2011.10.19 VIC Tamura */
SET_PROGRESS(100*(LOOP_COUNTER++)/max_size, 3, 1); // critical, but works
if(IS_CANSELED()) continue;

		const vector< VCData >::const_iterator it=m_vecdat.begin()+i_;
		const Double it_error = sqrt( it->Variance() * cv2 );
		const vector<VCData>::const_iterator it2_begin=lower_bound(m_vecdat.begin(), m_vecdat.end(), VCData::putVCData(it->Value()-it_error));
		const vector<VCData>::const_iterator it2_end=upper_bound(it2_begin, m_vecdat.end(), VCData::putVCData(it->Value()+it_error));

		for(vector<VCData>::const_iterator it2 = it2_begin; it2<it2_end; it2++)
		{
			const VCData diff = *it - *it2;
			if( diff.Value()*diff.Value() > cv2 * it2->Variance() ) continue;

			if( it2->putVecCoef().begin()->second == i3 )
			{
				budex.setIndex(-1, it2->putVecCoef().begin()->first,
								it->putVecCoef().begin()->first, it->putVecCoef().rbegin()->first);
			}
			else
			{
				budex.setIndex(-1, it2->putVecCoef().rbegin()->first,
								it->putVecCoef().begin()->first, it->putVecCoef().rbegin()->first);
			}
			if( budex.cross_product_412() <= 0.0 ) continue;

			if( it->putVecCoef().begin()->second == i3 )
			{
				budex2.setIndex(-1, it->putVecCoef().begin()->first,
								it2->putVecCoef().begin()->first, it2->putVecCoef().rbegin()->first);
			}
			else
			{
				budex2.setIndex(-1, it->putVecCoef().rbegin()->first,
								it2->putVecCoef().begin()->first, it2->putVecCoef().rbegin()->first);
			}

			if( budex2.cross_product_412() <= 0.0 ) continue;

#ifdef _OPENMP
				#pragma omp critical
#endif
			{
				if( found_bud_tray.find(budex) == found_bud_tray.end() )
				{
					ans.push_back( budex );
					found_bud_tray.insert( budex );
				}
				if( found_bud_tray.find(budex2) == found_bud_tray.end() )
				{
					ans.push_back( budex2 );
					found_bud_tray.insert( budex2 );
				}
			}
		}
	}
	}
	catch(bad_alloc& ball)
	{
		throw nerror(ball, __FILE__, __LINE__, __FUNCTION__);
	}

/*   2011.10.19 VIC Tamura */
CHECK_INTERRUPTION();
}

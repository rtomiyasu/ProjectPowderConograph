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
#include <list>
#include "indexing_func_dim2.hh"
#include "utility_data_structure/VCData.hh"
#include "utility_data_structure/Bud.hh"
#include "utility_data_structure/TreeLattice.hh"
#include "utility_func/stopx.hh"
#include "zlog/zlog.hh"
#include "GenerateRelation.hh"

bool find_bud(const vector<Bud>::const_iterator it_begin, 
		const vector<Bud>::const_iterator it_end, const Double& cv2,
		const Bud& budex)
{
	const Double Q1 = budex.Q1().Value();
	const Double Q1_var = cv2 * budex.Q1().Variance();
	
	const Double diff = sqrt(Q1_var);
	const vector<QData>& qdata = VCData::putPeakQData();
	
	const Int4 ibegin = distance( qdata.begin(), lower_bound( qdata.begin(), qdata.end(), QData(Q1 - diff, 0.0) ) );
	const Int4 iend = distance( qdata.begin(), upper_bound( qdata.begin()+ibegin, qdata.end(), QData(Q1 + diff, 0.0) ) );
	
	vector<Bud>::const_iterator it_basket;
	Bud budex2;
	for(Int4 i=ibegin; i<iend; i++)
	{
		budex2.setIndex( i, budex.iK2(), budex.iK3(), budex.iK4() );
		it_basket = lower_bound( it_begin, it_end, budex2 );
		
		if( it_basket < it_end && *it_basket == budex2 ) return true;
	}
	return false;
}

void findRelation(const Double& cv2, const Int4& max_peak_num, const Int4& max_itnum,
		vector<Bud>& ans)
{
	GenerateRelation gen_rel(max_peak_num);

	gen_rel.initialize_type2();
	gen_rel.putRelationCandidate_type2(cv2, ans);
	sort( ans.begin(), ans.end() );

ZLOG_INFO( "Number of (q1, q2, q3, q4) satisfying 2(q1 + q2) = q3 + q4 : " + num2str(ans.size()) + "\n" );
	
	const Int4 num_bud = ans.size();

	vector<Bud> bud_tray;
	gen_rel.initialize_type3();
	gen_rel.putRelationCandidate_type3(cv2, bud_tray);
	if( !bud_tray.empty() )
	{
		sort( bud_tray.begin(), bud_tray.end() );

		if( !find_bud(ans.begin(), ans.begin()+num_bud, cv2, bud_tray[0]) )
		{
			ans.push_back(bud_tray[0]);
		}

		const Int4 max_size = bud_tray.size();

#ifdef _OPENMP
		#pragma omp parallel for
#endif
		for(Int4 i_=1; i_<max_size; i_++)
		{
			const vector<Bud>::const_iterator it=bud_tray.begin()+i_;

			if( *(it-1) == *it ) continue;
			if( !find_bud(ans.begin(), ans.begin()+num_bud, cv2, *it) )
			{
#ifdef _OPENMP
				#pragma omp critical
#endif
				{
					ans.push_back(*it);
				}
			}
		}
	}

ZLOG_INFO( "Number of (q?, q2, q3, q4) detected from the equation 3q2 + q4 = 3q3 + q5 : " + num2str( ans.size() - num_bud ) + "\n" );

	sort( ans.begin(), ans.end() );
}




// In this function, budex.iK1() < 0 is assumed.
// On input and output, bud_basket is sorted.
static void find_new_root_unobs(const Double& cv2, const Bud budex,
		const vector<Bud>& bud_basket,
		set<Bud>& newly_added_peak)
{
	const VCData K4 = budex.Q1();
	const VCData K3 = ( budex.Q2() + budex.Q3() )* 2 - K4;

	const Double Q3 = K3.Value();
	const Double Q3_diff = sqrt( cv2 * K3.Variance() );
	
	const vector<QData>& qdata = VCData::putPeakQData();
	const pair<vector<QData>::const_iterator, vector<QData>::const_iterator>
		it_pair_dbl2 = closest_qdata(qdata.begin(), qdata.end(), Q3, Q3_diff);
	if( it_pair_dbl2.first >= it_pair_dbl2.second ) return;
	assert( qdata.begin() <= it_pair_dbl2.first );
	assert( it_pair_dbl2.second <= qdata.end() );

	const Int4 Q3_start = distance(qdata.begin(), it_pair_dbl2.first);
	const Int4 Q3_end = distance(qdata.begin(), it_pair_dbl2.second);
	
	Bud bud_tray = putBud124( budex.iK3(), budex.iK2(), -1 );

	pair<vector<Bud>::const_iterator, vector<Bud>::const_iterator> it_pair;
	vector<Bud>::const_iterator it2;

	Vec_BOOL found_flag(Q3_end-Q3_start, false);
	if( K3.Value() <= K4.Value() )
	{
		it_pair = equal_range_bud124(bud_basket.begin(), bud_basket.end(), bud_tray.iK1(), bud_tray.iK2(), bud_tray.iK4());
		for(it2=it_pair.first; it2<it_pair.second; it2++)
		{
			if( Q3_start <= it2->iK3() && it2->iK3() < Q3_end )
			{
				found_flag[it2->iK3() - Q3_start] = true;
			}
		}

		for(Int4 k2=Q3_start, k=0; k2<Q3_end; k2++, k++)
		{
			if( found_flag[k] ) continue;
			bud_tray.setIndex( bud_tray.iK1(), bud_tray.iK2(), k2, -1 );
//			bud_basket.insert( upper_bound( bud_basket.begin(), bud_basket.end(), bud_tray ), bud_tray );
			newly_added_peak.insert( bud_tray );
		}
	}
	else
	{
		it_pair = equal_range_bud12(bud_basket.begin(), bud_basket.end(), bud_tray.iK1(), bud_tray.iK2());

		vector<Bud>::const_iterator it, it2;
		for(it2=it_pair.first; it2<it_pair.second; it2++)
		{
			if( it2->iK3() == -1 && Q3_start <= it2->iK4() && it2->iK4() < Q3_end )
			{
				found_flag[it2->iK4() - Q3_start] = true;
			}
		}

		for(Int4 k2=Q3_start, k=0; k2<Q3_end; k2++, k++)
		{
			if( found_flag[k] ) continue;
			bud_tray.setIndex( bud_tray.iK1(), bud_tray.iK2(), -1, k2 );
//			bud_basket.insert( upper_bound( bud_basket.begin(), bud_basket.end(), bud_tray ), bud_tray );
			newly_added_peak.insert( bud_tray );
		}
	}
}


// In this function, K1, K2, K34 >= 0.
// On input and output, bud_basket is sorted.
static void find_new_root_obs(const Double& cv2, 
		const Int4& K1, const Int4& K2, const Int4& K34,
		const vector<Bud>& bud_basket,
		set<Bud>& newly_added_peak)
{
	const vector<QData>& qdata = VCData::putPeakQData();
	const Double Q3 = ( qdata[K1].q + qdata[K2].q ) * 2.0 - qdata[K34].q;
	if( Q3 > qdata[K34].q ) return;

	Bud bud_tray = putBud124( K1, K2, K34 );
	pair< vector<Bud>::const_iterator, vector<Bud>::const_iterator> it_pair = equal_range_bud124(bud_basket.begin(), bud_basket.end(), bud_tray.iK1(), bud_tray.iK2(), bud_tray.iK4());
	if( it_pair.first < it_pair.second ) return;

	newly_added_peak.insert( bud_tray );
}


// On input and output, bud_basket is sorted.
static void find_new_root(const Double& cv2, const Bud budex,
		const vector<Bud>& bud_basket,
		set<Bud>& newly_added_peak)
{
	if( budex.iK3() < 0 )
	{
		return;
	}

	if( budex.iK1() < 0 )
	{
		find_new_root_unobs(cv2, budex, bud_basket, newly_added_peak);
	}
	else // budex.iK3(), budex.iK1(), budex.iK2() >=0. 
	{
		find_new_root_obs(cv2, budex.iK3(), budex.iK1(), budex.iK2(), bud_basket, newly_added_peak);
		find_new_root_obs(cv2, budex.iK3(), budex.iK2(), budex.iK1(), bud_basket, newly_added_peak);
	}
}




// On input and output, bud_basket is sorted.
static void find_leftandright_branch(const Double& cv2,
		const vector<Bud>& bud_basket,
		Vec_BOOL& root_flag,
		vector< vector<Int4> >& leftbr_tray,
		vector< vector<Int4> >& rightbr_tray)
{
	root_flag.clear();
	leftbr_tray.clear();
	rightbr_tray.clear();
	
	const Int4 nsize = bud_basket.size();
	root_flag.resize(nsize, true);
	leftbr_tray.resize(nsize);
	rightbr_tray.resize(nsize);

	const vector<QData>& qdata = VCData::putPeakQData();
	const Int4 bad_num = bud_basket.size();

	Int4 k2;
	Double diff, Q3;
	Bud budex;
	pair<vector<Bud>::const_iterator, vector<Bud>::const_iterator> it_pair;

/*   2011.10.19 VIC Tamura */
Int4 LOOP_COUNTER = 0;

#ifdef _OPENMP
		#pragma omp parallel for private(k2, diff, it_pair, budex, Q3)
#endif
	for(Int4 k=0; k<bad_num; k++)
	{
/*   2011.10.19 VIC Tamura */
SET_PROGRESS(100*(LOOP_COUNTER++)/bad_num, 5, 1); // critical, but works
if(IS_CANSELED()) continue;
		
		const vector<Bud>::const_iterator it=bud_basket.begin() + k;
		if( it->iK3() == it->iK4() ) continue;

		if( it->Q2() < it->Q1() )
		{
			budex = putBud124(it->iK3(), it->iK2(), it->iK1());
			Q3 = ( it->Q2().Value() + it->Q3().Value() )* 2.0 - it->Q1().Value();
		}
		else
		{
			budex = putBud124(it->iK3(), it->iK1(), it->iK2());
			Q3 = ( it->Q1().Value() + it->Q3().Value() )* 2.0 - it->Q2().Value();
		}

		it_pair = equal_range_bud124(bud_basket.begin(), bud_basket.end(), budex.iK1(), budex.iK2(), budex.iK4());

		k2 = distance( bud_basket.begin(), it_pair.first );
		for(vector<Bud>::const_iterator it2=it_pair.first; it2<it_pair.second; it2++, k2++)
		{
			if( budex.iK1() < 0 || budex.iK2() < 0 || budex.iK4() < 0 )
			{
				diff = qdata[it2->iK3()].q - Q3;
				if( diff * diff > cv2 * qdata[it2->iK3()].q_var ) continue;
			}
			if( it->iK3() == budex.iK2() )
			{
#ifdef _OPENMP
				#pragma omp critical (lhs)
#endif
				{
					leftbr_tray[k2].push_back(k);
				}
			}
			if( it->iK3() == budex.iK1() )
			{
#ifdef _OPENMP
				#pragma omp critical (rhs)
#endif
				{
					rightbr_tray[k2].push_back(k);
				}
			}
			root_flag[k] = false;
		}
	}
/*   2011.10.19 VIC Tamura */
CHECK_INTERRUPTION();
}



// On input and output, bud_basket is sorted.
static void find_leftandright_root(const Double& cv2,
		const Bud& budex, const vector<Bud>& bud_basket,
		vector<Int4>& root_leftbr,
		vector<Int4>& root_rightbr)
{
	root_leftbr.clear();
	root_rightbr.clear();
	
	const vector<QData>& qdata = VCData::putPeakQData();
	
	Int4 k2;
	Double Q4, diff;
	vector<Bud>::const_iterator it2;
	pair<vector<Bud>::const_iterator, vector<Bud>::const_iterator> it2_pair;

	// Find a bud taking the root on the left branch.
	Bud bud_tray = putBud123( budex.iK2(), budex.iK3(), budex.iK1() );

	Q4 = ( budex.Q2().Value() + budex.Q3().Value() )*2.0 - budex.Q1().Value();
	if( Q4 > budex.Q1().Value() )
	{
		it2_pair = equal_range_bud12( bud_basket.begin(), bud_basket.end(), bud_tray.iK1(), bud_tray.iK2() );
		k2 = distance( bud_basket.begin(), it2_pair.first );
		for(it2=it2_pair.first; it2<it2_pair.second; it2++, k2++)
		{
			if( it2->iK3() != bud_tray.iK3() ) continue;
			if( it2->iK3() == it2->iK4() ) continue;
			if( bud_tray.iK1() < 0 || bud_tray.iK2() < 0 || bud_tray.iK3() < 0 )
			{
				diff = qdata[it2->iK4()].q - Q4;
				if( diff * diff > cv2 * qdata[it2->iK4()].q_var ) continue;
			}
			root_leftbr.push_back(k2);
		}
	}
	


	// Find a bud taking the root on the right branch.
	if( budex.iK1() == budex.iK2() )
	{
		root_rightbr = root_leftbr;
	}
	else
	{
		bud_tray = putBud123( budex.iK1(), budex.iK3(), budex.iK2() );

		Q4 = ( budex.Q1().Value() + budex.Q3().Value() )*2.0 - budex.Q2().Value();
		if( Q4 > budex.Q2().Value() )
		{
			it2_pair = equal_range_bud12( bud_basket.begin(), bud_basket.end(), bud_tray.iK1(), bud_tray.iK2() );
			k2 = distance( bud_basket.begin(), it2_pair.first );
			for(it2=it2_pair.first; it2<it2_pair.second; it2++, k2++)
			{
				if( it2->iK3() != bud_tray.iK3() ) continue;
				if( it2->iK3() == it2->iK4() ) continue;
				if( bud_tray.iK1() < 0 || bud_tray.iK2() < 0 || bud_tray.iK3() < 0 )
				{
//					if( it2->iK3() < 0 ) continue;
					diff = qdata[it2->iK4()].q - Q4;
					if( diff * diff > cv2 * qdata[it2->iK4()].q_var ) continue;
				}
				root_rightbr.push_back(k2);
			}
		}
	}
}


// On input, bud_basket is sorted.
void findNeighborRelation(const Double& cv2,
		vector<Bud>& bud_basket,
		Vec_BOOL& root_flag, 
		vector< vector<Int4> >& leftbr_tray,
		vector< vector<Int4> >& rightbr_tray,
		vector< vector<Int4> >& root_leftbr,
		vector< vector<Int4> >& root_rightbr)
{
	vector<Bud> new_bud_tray = bud_basket;
	set<Bud> new_bud_tray_all;

	do
	{
/*   2011.10.19 VIC Tamura */
SET_PROGRESS(100*(bud_basket.size()-new_bud_tray.size())/bud_basket.size(), 4, 1);
CHECK_INTERRUPTION();

		const Int4 isize = new_bud_tray.size();
#ifdef _OPENMP
		#pragma omp parallel
#endif
		{
			set<Bud> new_bud;

#ifdef _OPENMP
			#pragma omp for
#endif
			for(Int4 i=0; i<isize; i++)
			{
				find_new_root(cv2, new_bud_tray[i], bud_basket, new_bud);
			}

#ifdef _OPENMP
			#pragma omp single
#endif
			{
				new_bud_tray.clear();
			}

#ifdef _OPENMP
			#pragma omp critical
#endif
			{
				for(set<Bud>::const_iterator it=new_bud.begin(); it!=new_bud.end(); it++)
				{
					if( new_bud_tray_all.find(*it) != new_bud_tray_all.end() ) continue;
					new_bud_tray_all.insert(*it);
					new_bud_tray.push_back(*it);
				}
			}
		}
	} while( !new_bud_tray.empty() );

	bud_basket.insert(bud_basket.end(), new_bud_tray_all.begin(), new_bud_tray_all.end());
	sort( bud_basket.begin(), bud_basket.end() );

	const Int4 nsize = bud_basket.size();
	
	find_leftandright_branch(cv2, bud_basket, 
								root_flag, leftbr_tray, rightbr_tray);
	
	root_leftbr.clear();
	root_rightbr.clear();

	root_leftbr.resize(nsize);
	root_rightbr.resize(nsize);

/*   2011.10.19 VIC Tamura */
Int4 LOOP_COUNTER = 0;

#ifdef _OPENMP
	#pragma omp parallel for
#endif
	for(Int4 i=0; i<nsize; i++)
	{
/*   2011.10.19 VIC Tamura */
SET_PROGRESS(100*(LOOP_COUNTER++)/nsize, 6, 1); // critical, but works
if(IS_CANSELED()) continue;

		if( !root_flag[i] ) continue;
		if( bud_basket[i].iK3() == bud_basket[i].iK4() ) continue;

		find_leftandright_root(cv2, bud_basket[i], bud_basket, root_leftbr[i], root_rightbr[i]);
	}

/*   2011.10.19 VIC Tamura */
CHECK_INTERRUPTION();
}






static void construct_tree(const Int4& k, 
const vector<Bud>::const_iterator& it_begin,
const vector< vector<Int4> >& leftbr_tray,
const vector< vector<Int4> >& rightbr_tray,
TreeLattice& tree_tray)
{
	const vector<Bud>::const_iterator& it = it_begin + k;
	
	TreeLattice tree_left, tree_right;
	tree_left.setRoot(it->iK1(), it->iK4());
	tree_right.setRoot(it->iK2(), it->iK4());

	if( !leftbr_tray[k].empty() )
	{
		TreeLattice tree_left2;
		for(vector<Int4>::const_iterator it2=leftbr_tray[k].begin(); it2!=leftbr_tray[k].end(); it2++)
		{
			if( *it2 != k )
			{
				construct_tree(*it2, it_begin, leftbr_tray, rightbr_tray, tree_left2);
				if( tree_left.Root().count_bud() < tree_left2.Root().count_bud() )
				{
					tree_left = tree_left2;
				}
			}
		}
		if( tree_left.Root().Right() != it->iK4() ) tree_left.swapBranch();
	}
	
	if( !rightbr_tray[k].empty() )
	{
		TreeLattice tree_right2;
		for(vector<Int4>::const_iterator it2=rightbr_tray[k].begin(); it2!=rightbr_tray[k].end(); it2++)
		{
			if( *it2 != k )
			{
				construct_tree(*it2, it_begin, leftbr_tray, rightbr_tray, tree_right2);
				if( tree_right.Root().count_bud() < tree_right2.Root().count_bud() )
				{
					tree_right = tree_right2;
				}
			}
		}
		if( tree_right.Root().Right() != it->iK4() ) tree_right.swapBranch();
	}

	tree_tray.clear();
	tree_tray.setRoot(tree_left.Root(), tree_right.Root());
}


void constructTree(const vector<Bud>& bud_basket, 
const Vec_BOOL& root_flag, 
const vector< vector<Int4> >& leftbr_tray,
const vector< vector<Int4> >& rightbr_tray,
const vector< vector<Int4> >& root_leftbr_tray,
const vector< vector<Int4> >& root_rightbr_tray, const Int4& MAX_TREE,
vector<TreeLattice>& tree_tray)
{
	tree_tray.clear();
	const Int4 nsize = bud_basket.size();
	
/*   2011.10.19 VIC Tamura */
Int4 LOOP_COUNTER = 0;

#ifdef _OPENMP
	#pragma omp parallel
#endif
	{
		TreeLattice chibi_tree;
		TreeLattice tree_left, tree_left2, tree_right, tree_right2;
		Int4 index;
		multiset<TreeLattice> tree_tray2;

#ifdef _OPENMP
	#pragma omp for
#endif
		for(Int4 k=0; k<nsize; k++)
		{
/*   2011.10.19 VIC Tamura */
SET_PROGRESS(100*(LOOP_COUNTER++)/nsize, 7, 1); // critical, but works
if(IS_CANSELED()) continue;

			if( !root_flag[k] ) continue;
			const vector<Bud>::const_iterator it = bud_basket.begin() + k;

			const Int4 lsize = root_leftbr_tray[k].size();
			const Int4 rsize = root_rightbr_tray[k].size();

			// In the case of super-basis, avoid duplication of trees.
			for(index=0; index<lsize; index++)
			{
				if( k > root_leftbr_tray[k][index] ) break;
			}
			if( index<lsize ) continue;
			for(index=0; index<rsize;



					index++)
			{
				if( k > root_rightbr_tray[k][index] ) break;
			}
			if( index<rsize ) continue;
			construct_tree(k, bud_basket.begin(), leftbr_tray, rightbr_tray, chibi_tree);

			if( bud_basket[k].iK3() == bud_basket[k].iK4() )
			{
				chibi_tree.setRootEqualUpper();
			}
			else if( it->IsSuperBasisObtuse() )
			{
				tree_left.setRoot( it->iK2(), it->iK3() );
				for(vector<Int4>::const_iterator it2=root_leftbr_tray[k].begin(); it2<root_leftbr_tray[k].end(); it2++)
				{
					construct_tree(*it2, bud_basket.begin(), leftbr_tray, rightbr_tray, tree_left2);
					if( tree_left.Root().count_bud() < tree_left2.Root().count_bud() )
					{
						tree_left = tree_left2;
					}
				}
				chibi_tree.setRootOnLeftBranch(tree_left.Root());

				tree_right.setRoot( it->iK3(), it->iK1() );
				for(vector<Int4>::const_iterator it2=root_rightbr_tray[k].begin(); it2<root_rightbr_tray[k].end(); it2++)
				{
					construct_tree(*it2, bud_basket.begin(), leftbr_tray, rightbr_tray, tree_right2);
					if( tree_right.Root().count_bud() < tree_right2.Root().count_bud() )
					{
						tree_right = tree_right2;
					}
				}
				chibi_tree.setRootOnRightBranch(tree_right.Root());
			}
			else if( it->Q1() < it->Q2() )
			{
				NodeB node = chibi_tree.Root();
				if( lsize != 0 )
				{
					construct_tree(root_leftbr_tray[k][0], bud_basket.begin(), leftbr_tray, rightbr_tray, tree_left);
					for(vector<Int4>::const_iterator it2=root_leftbr_tray[k].begin()+1; it2<root_leftbr_tray[k].end(); it2++)
					{
						construct_tree(*it2, bud_basket.begin(), leftbr_tray, rightbr_tray, tree_left2);
						if( tree_left.Root().count_bud() < tree_left2.Root().count_bud() )
						{
							tree_left = tree_left2;
						}
					}
					if( tree_left.Root().Right() != it->iK2() ) tree_left.swapBranch();
					chibi_tree.setRoot(node, tree_left.Root());
				}
				else
				{
					const NodeB node2(it->iK3(), it->iK2());
					chibi_tree.setRoot(node, node2);
				}
			}
			else
			{
				NodeB nodex = chibi_tree.Root();
				nodex.swapBranch();
				if( rsize != 0 )
				{
					construct_tree(root_rightbr_tray[k][0], bud_basket.begin(), leftbr_tray, rightbr_tray, tree_right);
					for(vector<Int4>::const_iterator it2=root_rightbr_tray[k].begin()+1; it2<root_rightbr_tray[k].end(); it2++)
					{
						construct_tree(*it2, bud_basket.begin(), leftbr_tray, rightbr_tray, tree_right2);
						if( tree_right.Root().count_bud() < tree_right2.Root().count_bud() )
						{
							tree_right = tree_right2;
						}
					}
					if( tree_right.Root().Right() != it->iK1() ) tree_right.swapBranch();
					chibi_tree.setRoot(tree_right.Root(), nodex);
				}
				else
				{
					const NodeB nodex2(it->iK3(), it->iK1());
					chibi_tree.setRoot(nodex2, nodex);
				}
			}

			chibi_tree.setSortingCriteria();
			tree_tray2.insert(chibi_tree);
			if( (Int4)tree_tray2.size() > MAX_TREE ) tree_tray2.erase(--tree_tray2.end());
		}
#ifdef _OPENMP
		#pragma omp critical
#endif
		{
			tree_tray.insert(tree_tray.end(), tree_tray2.begin(), tree_tray2.end());
		}
	}
/*   2011.10.19 VIC Tamura */
CHECK_INTERRUPTION();

	sort(tree_tray.begin(), tree_tray.end());
	if( (Int4)tree_tray.size() > MAX_TREE ) tree_tray.resize(MAX_TREE);
}

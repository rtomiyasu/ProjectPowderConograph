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
#include <time.h>
#include "utility_data_structure/index_set.hh"
#include "utility_data_structure/Bud.hh"
#include "utility_data_structure/Bud2.hh"
#include "utility_data_structure/Node3.hh"
#include "utility_data_structure/TreeLattice.hh"
#include "utility_lattice_reduction/put_Selling_reduced_lattice.hh"
#include "PeakPosData.hh"
#include "ControlParam.hh"
#include "p_out_indexing.hh"
#include "zerror_type/error_out.hh"
#include "indexing_func_dim2.hh"
#include "indexing_func_dim3.hh"
#include "IndexingLattice.hh"
#include "model_function/PeakPosModel.hh"
#include "zlog/zlog.hh"
#include "zparam/ZParawError.hh"
#include "utility_func/stopx.hh"
#include "lattice_symmetry/VCLatticeFigureOfMeritToCheckSymmetry.hh"


IndexingLattice::IndexingLattice()
{
   	m_max_peak_num = 30;
	m_max_edge_num = 300;
	m_max_node_num = 300;

	m_min_unitcell_volume = 5.0;
	m_max_unitcell_volume = 10000.0;
   	m_min_primitive_cell_edge = 2.0;
   	m_cv2 = 0.0;

	m_etype_peak_shift = kPeakShiftFunction_Type0;
	m_WlengthX = 1.54056;

	m_search_level = ConographQuickSearch;
}


IndexingLattice::~IndexingLattice()
{
}


// Set the member variables.
void IndexingLattice::setParam(const ControlParam& cont) 
{
   	m_max_peak_num = cont.putMaxPeakNum();
   	m_max_edge_num = cont.putMaxEdgeNum();
   	m_max_node_num = cont.putMaxNodeNum();

   	m_min_unitcell_volume = cont.putMinUnitCellVolume();
   	m_max_unitcell_volume = cont.putMaxUnitCellVolume();
   	m_min_primitive_cell_edge = cont.putMinLatticePointDistance();
 	m_cv2 = cont.putCriticalValueSquareForLinearSum();

	m_etype_peak_shift = cont.putPeakShiftFunctionType();
	m_WlengthX = cont.putWaveLength();

	const vector<Double>& peak_shift_param_rad = cont.putPeakShiftParamRadian();
	const Int4 param_num = peak_shift_param_rad.size();
	assert( m_etype_peak_shift != kPeakShiftFunction_Type0 || param_num == 0 );
	assert( m_etype_peak_shift != kPeakShiftFunction_Type1 || param_num == 1 );

	m_peak_shift_param_rad.resize(param_num);
	for(Int4 i=0; i<param_num; i++)
	{
		m_peak_shift_param_rad[i] = peak_shift_param_rad[i];
	}

	m_search_level = cont.putSearchLevel();
}





static void addBudBranches(const Int4& index,
		const vector< vector<Int4> >& leftbr_tray, const vector< vector<Int4> >& rightbr_tray,
		set<Int4>& selected_bud_index)
{
	selected_bud_index.insert(index);

	for(UInt4 l=0; l<leftbr_tray.at(index).size(); l++)
	{
		addBudBranches(leftbr_tray[index][l], leftbr_tray, rightbr_tray, selected_bud_index);
	}

	for(UInt4 l=0; l<rightbr_tray.at(index).size(); l++)
	{
		addBudBranches(rightbr_tray[index][l], leftbr_tray, rightbr_tray, selected_bud_index);
	}
}



void IndexingLattice::determineTwoDimLattices(
		const PeakPosData& pdata,
		const string& fname)
{
	static const Int4 MAX_ITNUM = 100000000;
	
ZLOG_INFO( "\nSearching for q-values satisfying the Ito equation 2(q1 + q2) = q3 + q4...\n" );
	vector<Bud> bud_basket;
	findRelation(m_cv2, m_max_peak_num, MAX_ITNUM, bud_basket);

if( bud_basket.empty() )
{
	throw ZErrorMessage("No q-values satisfy the Ito equation.  Stop.\n", __FILE__, __LINE__, __FUNCTION__);
}

	Vec_BOOL root_flag2;
	vector< vector<Int4> > leftbr_tray, rightbr_tray, root_leftbr, root_rightbr;

	const Int4 num_bud_ini = bud_basket.size();
	findNeighborRelation(m_cv2, bud_basket, root_flag2,
									leftbr_tray, rightbr_tray, root_leftbr, root_rightbr);
ZLOG_INFO( "Number of (q1, q2, q3, q?) added under supposition of extinction rules: " + num2str( bud_basket.size() - num_bud_ini ) + "\n\n" );

#ifdef DEBUG
ZLOG_INFO( "Outputting a Bud-file...\n" );
	print_bud_data(bud_basket, root_flag2,
		 				leftbr_tray, rightbr_tray, root_leftbr, root_rightbr, fname+".bud");
#endif

//	static const Int4 num_topograph_to_fit_zero_shift = 10;
	if( (Int4)bud_basket.size() <= m_max_edge_num )
	{
		if( m_etype_peak_shift == kPeakShiftFunction_Type0 )
		{
			m_selected_bud_tray = bud_basket;
			return;
		}
	}
//	const Int4 max_edge_num = ((Int4)bud_basket.size()>m_max_edge_num?m_max_edge_num:num_topograph_to_fit_zero_shift);

ZLOG_INFO( "Making 2-dimensional topographs...\n" );
	vector<TreeLattice> tree_tray;
	constructTree(bud_basket, root_flag2,
					leftbr_tray, rightbr_tray, root_leftbr, root_rightbr, 
					m_max_edge_num, tree_tray);

ZLOG_INFO( "Number of topographs constructed : " + num2str( tree_tray.size() ) + "\n" );

#ifdef DEBUG
ZLOG_INFO( "Outputting a Tree-file...\n" );
	TreeLattice::print(fname+".tree", tree_tray.begin(), tree_tray.end() );
#endif

	m_selected_bud_tray.clear();
	if( tree_tray.empty() ) return;

	set<Int4> selected_bud_index;
	set<Bud> tray2;
	vector<Bud>::const_iterator it_basket;

	Int4 last_tree = tree_tray.size()-1;
	for(vector<TreeLattice>::const_iterator it=tree_tray.begin(); it<tree_tray.end(); it++)
	{
		it->putRootBuds(tray2);
		for(set<Bud>::const_iterator it2=tray2.begin(); it2!=tray2.end(); it2++)
		{
			it_basket = lower_bound(bud_basket.begin(), bud_basket.end(), *it2);
			if( it_basket >= bud_basket.end() || !(*it_basket == *it2) )
			{
				continue;
			}
			assert( it_basket + 1 >= bud_basket.end() || !( *(it_basket+1) == *it2 ) );

			addBudBranches( distance((vector<Bud>::const_iterator)bud_basket.begin(), it_basket),
							leftbr_tray, rightbr_tray, selected_bud_index);
		}

		if( (Int4)selected_bud_index.size() >= m_max_edge_num )
		{
			last_tree = distance((vector<TreeLattice>::const_iterator)tree_tray.begin(), it)+1;
			break;
		}
	}

ZLOG_INFO( "Maximum number of q_values contained in a topograph : " + num2str( tree_tray.begin()->putCountOfQ() ) + "\n"
			+ "Number of q_values contained in the " + num2str<Int4>( last_tree + 1 ) + "th topograph : " + num2str( (tree_tray.begin() + last_tree)->putCountOfQ() ) + "\n"
			+ "The program selected " + num2str( selected_bud_index.size() )
			+ " quadruples (q1, q2, q3, q4) from " + num2str<Int4>( last_tree + 1 ) + " topographs.\n\n" );

	Int4 index = 0;
	for(set<Int4>::const_iterator it=selected_bud_index.begin(); it!=selected_bud_index.end(); it++, index++)
	{
		m_selected_bud_tray.push_back(bud_basket[*it]);
	}

#ifdef DEBUG
ZLOG_INFO( "Outputting a Bud-file...\n" );
	sort( m_selected_bud_tray.begin(), m_selected_bud_tray.end() );
	print_bud_data(m_selected_bud_tray, fname+".bud2");
#endif

}




template <class T>
static bool cmp_size(const index_set< vector<T> >& lhs, const index_set< vector<T> >& rhs)
{
	return lhs.element.size() < rhs.element.size();
}



// On input, bud_basket and node_basket ares sorted.
void IndexingLattice::determineThreeDimLattices(// const Double& MIN_NormM, const Double& MIN_RevM,
		const string& fname0)
{
	m_S_super.clear();
	if( m_selected_bud_tray.empty() ) return;

	set<Bud2> bud_tray;
	Bud2 nodex;
	for(vector<Bud>::const_iterator it2=m_selected_bud_tray.begin(); it2<m_selected_bud_tray.end(); it2++)
	{
		if( it2->iK3() < 0 || it2->iK4() < 0 ) continue;

		nodex.setIndex( it2->iK1(), it2->iK2(), it2->iK4(), it2->iK3() );
		bud_tray.insert( nodex );

		nodex.setIndex( it2->iK2(), it2->iK1(), it2->iK4(), it2->iK3() );
		bud_tray.insert( nodex );

		nodex.setIndex( it2->iK1(), it2->iK2(), it2->iK3(), it2->iK4() );
		bud_tray.insert( nodex );

		nodex.setIndex( it2->iK2(), it2->iK1(), it2->iK3(), it2->iK4() );
		bud_tray.insert( nodex );
	}

	vector<Bud2> bud_basket2( bud_tray.begin(), bud_tray.end() );

ZLOG_INFO( "The number of nodes (q1, q2, q3), (q1, q2, q4) generated from quadruples (q1, q2, q3, q4) : " + num2str<Int4>( bud_basket2.size() ) + "\n"
			+ "Searching for combinations of (q1, q2, q3), (q1, q4, q5) and q6 to construct 3-dim lattices....(This process might take several minutes...)\n" );
	clock_t start;
	start = clock();    /* Record computational time */

	const Int4 num_ref_figure_of_merit = min((Int4)m_num_ref_figure_of_merit, (Int4)VCData::putPeakQData().size());
	set< pair<Double, SymMat43_VCData> > node_basket_vol;
	set< pair<Double, SymMat43_VCData> > node_basket_FOM;
	combinate_bud(bud_basket2, m_cv2, m_min_unitcell_volume, m_max_unitcell_volume,
					m_min_primitive_cell_edge * m_min_primitive_cell_edge,
					num_ref_figure_of_merit, 1.0,
				    m_etype_peak_shift, m_WlengthX, m_peak_shift_param_rad,
				    (m_search_level==ConographQuickSearch?m_max_node_num:m_max_node_num/2),
				    m_search_level, node_basket_vol, node_basket_FOM);

	const Double calculation_time = (clock() - start) / CLOCKS_PER_SEC;
ZLOG_INFO( "Number of solutions generated from combinations of four nodes (q1,q2,q3) in CPU time "
			+ num2str( calculation_time ) + " [sec.] : " + num2str( node_basket_vol.size() + node_basket_FOM.size() ) + "\nSolutions satisfy one of the following conditions.\n" );
	if( !node_basket_vol.empty() )
	{
ZLOG_INFO( "> Unitcell volume : " + num2str( sqrt(1.0 / node_basket_vol.rbegin()->first) ) + "--"
				+ num2str( 1.0 / sqrt(node_basket_vol.begin()->first ) ) + "\n" );
	}
	if( !node_basket_FOM.empty() )
	{
ZLOG_INFO( "> Figure of merit (" + putLabel(SCWuM) + num2str(num_ref_figure_of_merit) + ") : " + num2str( -node_basket_FOM.rbegin()->first ) + "--"
				+ num2str( -node_basket_FOM.begin()->first ) + "\n" );
	}
ZLOG_INFO( "\n" );

	for(set< pair<Double, SymMat43_VCData> >::const_iterator it=node_basket_vol.begin(); it!=node_basket_vol.end(); it++)
	{
		m_S_super.push_back(it->second);
	}
	for(set< pair<Double, SymMat43_VCData> >::const_iterator it=node_basket_FOM.begin(); it!=node_basket_FOM.end(); it++)
	{
		m_S_super.push_back(it->second);
	}
}

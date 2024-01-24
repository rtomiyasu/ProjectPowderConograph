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
#include "../ControlParam.hh"
#include "OutputInfo.hh"

static vector<bool(*)(const VCLatticeFigureOfMeritToCheckSymmetry&, const VCLatticeFigureOfMeritToCheckSymmetry&)> putCmpFunctions()
{
	vector<bool(*)(const VCLatticeFigureOfMeritToCheckSymmetry&, const VCLatticeFigureOfMeritToCheckSymmetry&)> ans;
	ans.push_back( VCLatticeFigureOfMeritToCheckSymmetry::cmpFOMdeWolff );
	ans.push_back( VCLatticeFigureOfMeritToCheckSymmetry::cmpFOMWu );
	ans.push_back( VCLatticeFigureOfMeritToCheckSymmetry::cmpFOMReversed );
	ans.push_back( VCLatticeFigureOfMeritToCheckSymmetry::cmpFOMSymmetric );
	ans.push_back( VCLatticeFigureOfMeritToCheckSymmetry::cmpNumberOfNeighbors );

	assert(ans.size() == (size_t)putNumberOfSortCriterion());
	return ans;
}

const vector<bool(*)(const VCLatticeFigureOfMeritToCheckSymmetry&, const VCLatticeFigureOfMeritToCheckSymmetry&)> OutputInfo::CmpFunc = putCmpFunctions();

OutputInfo::OutputInfo()
{  
	candidate_num = 0;
	selected_candidate_label_all.resize(putNumberOfSortCriterion(), -1);
	selected_candidate_label_list.resize(putNumberOfSortCriterion(), -1);
};

void OutputInfo::setLabel(
		const vector<VCLatticeFigureOfMeritToCheckSymmetry>& lattice_result,
		const ControlParam& cData)
{
	this->lattice_index.clear();
	this->candidate_num = 0;

	const size_t lsize = lattice_result.size();
	const size_t SortCriterionNum = putNumberOfSortCriterion();

	vector<Int4> selected_candidate_index_all(SortCriterionNum, -1);
	vector<Int4> selected_candidate_index_list(SortCriterionNum, -1);
	VecDat3<Double> length_axis, angle_axis;

	for(size_t l=0; l<lsize; l++)
	{
		bool& output_flag = lattice_index.insert(pair<lattice_label, pair<Int4, bool> >(lattice_result[l].putLatticeLabel(), pair<Int4, bool>(l, false))).first->second.second;

		const LatticeFigureOfMerit& latfom = lattice_result[l].putLatticeFigureOfMerit();

		if( !cData.putOutputSymmetry( latfom.enumBravaisType() ) ) continue;
		if( latfom.putFiguresOfMerit().putContinuousNumberOfHKLInRange() > Double( cData.putMaxNumPeakInRange() ) ) continue;
		if( latfom.putFiguresOfMerit().putContinuousNumberOfHKLInRange() < Double( cData.putMinNumPeakInRange() ) ) continue;
		const Double& minABC = cData.putMinUnitCellEdgeABC();
		const Double& maxABC = cData.putMaxUnitCellEdgeABC();
		latfom.putReducedLatticeConstantsDegree(cData.putBaseCenteredAxis(), cData.putRhombohedralAxis(), length_axis, angle_axis);
		if( length_axis[0] < minABC || length_axis[1] < minABC || length_axis[2] < minABC ) continue;
		if( length_axis[0] > maxABC || length_axis[1] > maxABC || length_axis[2] > maxABC ) continue;

		if( latfom.putFiguresOfMerit().putFigureOfMeritWolff() < cData.putMinFOM() ) continue;

		for(size_t i=0; i<SortCriterionNum; i++)
		{
			if( selected_candidate_index_all[i] < 0
				|| OutputInfo::CmpFunc[i]( lattice_result[l], lattice_result[selected_candidate_index_all[i]] ) )
			{
				selected_candidate_index_all[i] = l;
			}
		}

		if( lattice_result[l].putNumberOfLatticesInNeighborhood() < 0 ) continue;

		output_flag = true;
		this->candidate_num++;

		for(size_t i=0; i<SortCriterionNum; i++)
		{
			if( selected_candidate_index_list[i] < 0
				|| OutputInfo::CmpFunc[i]( lattice_result[l], lattice_result[selected_candidate_index_list[i]] ) )
			{
				selected_candidate_index_list[i] = l;
			}
		}
	}

	for(size_t i=0; i<SortCriterionNum; i++)
	{
		if( selected_candidate_index_all[i] < 0 )
		{
			this->selected_candidate_label_all[i] = -1;
		}
		else
		{
			this->selected_candidate_label_all[i] = lattice_result[selected_candidate_index_all[i]].putLatticeLabel();
		}
		if( selected_candidate_index_list[i] < 0 )
		{
			this->selected_candidate_label_list[i] = -1;
		}
		else
		{
			this->selected_candidate_label_list[i] = lattice_result[selected_candidate_index_list[i]].putLatticeLabel();
		}
	}

	if( selected_candidate_index_all[(size_t)SCNN] >= 0
		&& lattice_result[selected_candidate_index_all[(size_t)SCNN]].putNumberOfLatticesInNeighborhood() < 1 )
	{
		this->selected_candidate_label_all[(size_t)SCNN] = -1;
	}
	if( selected_candidate_index_list[(size_t)SCNN] >= 0
			&& lattice_result[selected_candidate_index_list[(size_t)SCNN]].putNumberOfLatticesInNeighborhood() < 1 )
	{
		this->selected_candidate_label_list[(size_t)SCNN] = -1;
	}
}


pair<eBravaisType, lattice_label> put_selected_lattice_from_all(
			const vector<VCLatticeFigureOfMeritToCheckSymmetry> lattice_result[],
			const OutputInfo outinfo[], 
			const eSortCriterion& sort_criterion)
{
	bool (*cmp_func)(const VCLatticeFigureOfMeritToCheckSymmetry&, const VCLatticeFigureOfMeritToCheckSymmetry&) = OutputInfo::CmpFunc[(size_t)sort_criterion];

	pair<eBravaisType, lattice_label> ans;
	ans.first = Triclinic;
	ans.second = -1;
	const VCLatticeFigureOfMeritToCheckSymmetry* ans_lattice = NULL;
	
	static const Int4 NUM_LS = put_number_of_bravais_types();
	for(Int4 i=0; i<NUM_LS; i++)
	{
		const Int4 index = outinfo[i].putIndex( outinfo[i].putLatticeLabelSelectedAmongAll(sort_criterion) );
		if( index >= 0 ) 
		{
			if( ans.second < 0 || cmp_func( lattice_result[i][index], *ans_lattice ) )
			{
				ans.first = eBravaisType(i);
				ans.second = lattice_result[i][index].putLatticeLabel();
				ans_lattice = &lattice_result[i][index];
			}
		}
	}
	return ans;
}

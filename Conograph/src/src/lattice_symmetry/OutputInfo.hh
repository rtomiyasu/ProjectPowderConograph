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
#ifndef OUTPUTINFO_HH_
#define OUTPUTINFO_HH_

#include <map>
#include "../RietveldAnalysisTypes.hh"
#include "../bravais_type/enumBravaisType.hh"
#include "enumSortCriterion.hh"
#include "VCLatticeFigureOfMeritToCheckSymmetry.hh"

class ControlParam;

class OutputInfo
{
private:
	// The second entry indicates which index the lattice of the label corresponds to,
	// and whether the lattice should be output.
	map<lattice_label, pair<Int4, bool> > lattice_index;
	Int4 candidate_num; 					// Number of candidates to output for each crystal system.

	vector<lattice_label> selected_candidate_label_all;	// Solutions with the best score for each sorting criterion.
	vector<lattice_label> selected_candidate_label_list;	// Solutions with the best score contained in the list to output.
	
public:
	OutputInfo();

	inline Int4 putIndex(const lattice_label& arg) const;
	inline bool IsOutput(const lattice_label& arg) const;
	inline const lattice_label& putLatticeLabelSelectedAmongAll(const eSortCriterion& i) const { return selected_candidate_label_all[(size_t)i]; };
	inline const lattice_label& putLatticeLabelSelectedFromListToOutput(const eSortCriterion& i) const { return selected_candidate_label_list[(size_t)i]; };
	inline const Int4& putCandidateNumToOutput() const { return candidate_num; };

	void setLabel(const vector<VCLatticeFigureOfMeritToCheckSymmetry>& lattice_result,
					const ControlParam& cData);

	static const vector<bool(*)(const VCLatticeFigureOfMeritToCheckSymmetry&, const VCLatticeFigureOfMeritToCheckSymmetry&)> CmpFunc;
};

inline Int4 OutputInfo::putIndex(const lattice_label& arg) const
{
	map<lattice_label, pair<Int4, bool> >::const_iterator it=lattice_index.find(arg);
	if( it == lattice_index.end() ) return -1;
	else return it->second.first;
};

inline bool OutputInfo::IsOutput(const lattice_label& arg) const
{
	map<lattice_label, pair<Int4, bool> >::const_iterator it=lattice_index.find(arg);
	if( it == lattice_index.end() ) return false;
	else return it->second.second;
};

// Returns the smallest candidate and its crystal system 
// when outinfo[i].selected_lattice for all the crystal systems are compared by cmpfunc.
pair<eBravaisType, lattice_label> put_selected_lattice_from_all(
			const vector<VCLatticeFigureOfMeritToCheckSymmetry> lattice_result[],
			const OutputInfo outinfo[], 
			const eSortCriterion& sort_criterion);

#endif /*OUTPUTINFO_HH_*/

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
#ifndef coset_representative_data_HH_
#define coset_representative_data_HH_

#include<set>
#include"../RietveldAnalysisTypes.hh"
#include"enumPointGroup.hh"
#include"enumGroupToMaxSubgp.hh"
#include"../symmetric_operation/SymmetricOperation.hh"
#include"../symmetric_operation/XYZCoord.hh"
#include"../symmetric_operation/S1.hh"

const ePointGroup& enumLowerGroup(const eGroupToMaxSubgp&);
const ePointGroup& enumUpperGroup(const eGroupToMaxSubgp&);
const Int4& Index(const eGroupToMaxSubgp&);

eGroupToMaxSubgp enumMaxNormalSubgroup(const ePointGroup&);

void enumMaxSubgroup(const ePointGroup&, vector<eGroupToMaxSubgp>&);
void enumMinUpperGroup(const ePointGroup&, vector<eGroupToMaxSubgp>&);

// Returns the representatives of the equivalent classes which do not include the unit element. 
bool CosetRepresentativeMaxSubgp(const eGroupToMaxSubgp&, Int4&, eSymmetricOperation&);

// Returns true if and only if egp is a subgroup of egp2.
bool IsSubgroup(const ePointGroup& epg, const ePointGroup& epg2);

bool generateGroup(const ePointGroup&, const ePointGroup&, ePointGroup&);
bool generateGroup(const vector<ePointGroup>&, ePointGroup&);
bool generateGroup(const vector<SymmetricOperation>&, ePointGroup&);

// Returns all the subgroups of epg.
void enumSubgroup(const ePointGroup& epg, set<ePointGroup>&);

// Returns all the subgroups of epg that contains esubgp and does not equal esubgp.
void enumSubgroup(const ePointGroup& egp, const ePointGroup& esubgp, set<ePointGroup>&);

//void putGenerator(const ePointGroup& epg, vector<eSymmetricOperation>& soptray);

#endif /*coset_representative_data_HH_*/

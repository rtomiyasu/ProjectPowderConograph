/*
 * The MIT License

   LatticeWithSameQ (a module in Conograph to obtain all the lattices with the same computed lines)

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
#ifndef _LatticeWithSameQ_hh_
#define _LatticeWithSameQ_hh_
// LatticeWithSameQ.hh

#include "../RietveldAnalysisTypes.hh"

class BravaisType;
class ControlParam;
class LatticeFigureOfMerit;
class LatticeMetricTensor;
class PeakPosData;
class ReducedLatticeToCheckBravais;

using namespace std;


class LatticeWithSameQ
{
private:
	enum{ m_scale_of_superlattice = 1 };

	vector<Double> m_set_of_qvalues;
	Int4 m_dimension_of_lattice;
	Double m_resol;
	vector<LatticeMetricTensor> m_result;

	// If a unitcell in m_result is almost equivalent to arg, -1 is set to the unitcell.
	// Otherwise, if the unitcell has the best figure of merit among the neighbor unitcells,
	// the number of the neighbor unitcells is set.
	// To the other neighbor unitcells, -1 is set.
	// After this method, the entries of m_result are sorted by the unit-cell volume.
	void setNumberOfNeighbors(const SymMat<Double>& arg);

public:
    LatticeWithSameQ();
    ~LatticeWithSameQ();

    void setParam(const ControlParam&, const SymMat<Double>& S_super);

    // On output,
    // evaluation_ans.first: all the solusions were gained?
    // evaluation_ans.second: max of q-values used to determine the diagonal entries of solutions.
    // The second argument arg is used in the member method setNumberOfNeighbors.
    void execute(pair<bool, Double>& evaluation_ans, const SymMat<Double>& arg=0);

    // S_super is assumed to be Delone reduced.
    // Returns almost equivalent Bravais lattices after sorting them by the distance from S_super.
    inline const vector<Double>& putSetOfQ() const { return m_set_of_qvalues; };
    inline const vector<LatticeMetricTensor>& putResults() const { return m_result; };
    void putLatticesWithSameQ(const ControlParam& cdata,
    							const vector<QData>& set_of_qvalues,
    							vector<LatticeFigureOfMerit>& tray) const;
};

#endif

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
#ifndef _SortingLattice_hh_
#define _SortingLattice_hh_
// SortingLattice.hh

#include "RietveldAnalysisTypes.hh"
#include "utility_data_structure/SymMat43_VCData.hh"
#include "lattice_symmetry/ReducedVCLatticeToCheckBravais.hh"

class BravaisType;
class ControlParam;
class VCLatticeFigureOfMeritToCheckSymmetry;
class PeakPosData;

using namespace std;


class SortingLattice
{
private:
	enum{ NUM_LS = 14 };	// 	Triclinic, Monoclinic, Monoclinic(C),
							// Orthorhombic(P), Orthorhombic(C), Orthorhombic(I), Orthorhombic(F),
							// Tetragonal(P), Tetragonal(I), Rhombohedral, Hexiagonl
							// Cubic(P), Cubic(I), Cubic(F).

	static const bool m_DoesPrudentSymSearch;
	bool OutputSymmetry[NUM_LS];
	bool JudgeSymmetry[NUM_LS];

    Int4 m_num_ref_figure_of_merit;
    ePeakShiftFunctionType m_etype_peak_shift;
    Double m_WlengthX;
    vector<ZParawError> m_peak_shift_param_rad;

    static const Double m_cv2;
    Double m_resol;

	void putCentringTypes(const ReducedVCLatticeToCheckBravais& RLCB,
										const VCLatticeFigureOfMeritToCheckSymmetry& lattice_original,
										const BravaisType& brat,
										vector<VCLatticeFigureOfMeritToCheckSymmetry>& lattice_result) const;

	void putLatticeCandidatesForTriclinic(const vector<SymMat43_VCData>& S_super,
									const Double& MIN_NormM,
									const Double& MIN_RevM,
									vector<VCLatticeFigureOfMeritToCheckSymmetry>& result) const;

	void putLatticeCandidatesForEachBravaisTypes(
								const Double& MIN_NormM,
								const Double& MIN_RevM,
								const Int4& MAX_SIZE,
								const eABCaxis& baxis_flag,
								const eRHaxis& rhom_flag,
								vector<VCLatticeFigureOfMeritToCheckSymmetry> result[NUM_LS]) const;

public:
    SortingLattice();
    ~SortingLattice();

    void setParam(const ControlParam&);

	void putLatticeCandidatesForEachBravaisTypes(const vector<SymMat43_VCData>& S_super,
								const Double& MIN_NormM,
								const Double& MIN_RevM,
								const Int4& MAX_SIZE,
								const eABCaxis& baxis_flag,
								const eRHaxis& rhom_flag,
								vector<VCLatticeFigureOfMeritToCheckSymmetry> result[NUM_LS]) const;

	// The variables VCLatticeFigureOfMeritToCheckSymmetry::num_lattice_found in the argument are set.
	// For lattices having the best figure of merit (comparison by CmpFunc ) among solutions in the neighborhood,
	// the number of the solution in the neighborhood is set.
	// For other lattices, -1 is set. (This is for convenience in the class of OutputInfo.)
	// After this method, each result [i] is sorted by the unit-cell volume.
	void setNumberOfNeighbors(const eABCaxis& baxis_flag,
			bool (*CmpFunc)(const VCLatticeFigureOfMeritToCheckSymmetry&, const VCLatticeFigureOfMeritToCheckSymmetry&),
			vector<VCLatticeFigureOfMeritToCheckSymmetry> result[NUM_LS]) const;
};

#endif

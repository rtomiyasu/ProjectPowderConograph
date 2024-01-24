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
#ifndef VCLATTICEFIGUREOFMERITTOCHECKSYMMETRY_HH_
#define VCLATTICEFIGUREOFMERITTOCHECKSYMMETRY_HH_

#include "LatticeFigureOfMeritZeroShift.hh"
#include "../utility_data_structure/SymMat43_VCData.hh"

class OutputInfo;
typedef Int4 lattice_label;

// Class for outputting information about a lattice in index file.
class VCLatticeFigureOfMeritToCheckSymmetry
{
private:
	enum { NUM_LS = 14 };
	static const string CS_LABEL[NUM_LS];
	lattice_label m_label;
	
	LatticeFigureOfMeritZeroShift m_latfom;

	// m_S_red.first is Buerger-reduced.
	// m_S_red.second * m_S_red.first * Transpose(m_S_red.second) is obtuse.
	SymMat43_VCData m_S_red;
	Int4 m_num_lattice_found;

	bool checkIfLatticeIsMonoclinic(const ePointGroup& epg_new, const Double& cv2,
								map< SymMat<VCData>, NRMat<Int4> >& ans) const;
	bool checkIfLatticeIsOrthorhombic(const Double& cv2,
								map< SymMat<VCData>, NRMat<Int4> >& ans) const;
	bool checkIfLatticeIsTetragonal(const Double& cv2,
								map< SymMat<VCData>, NRMat<Int4> >& ans) const;
	bool checkIfLatticeIsHexagonal(const ePointGroup& epg_new, const Double& cv2,
								map< SymMat<VCData>, NRMat<Int4> >& ans) const;

	void setLatticeConstants43(const BravaisType& brat, const SymMat43_VCData& S);

public:
//	VCLatticeFigureOfMeritToCheckSymmetry();
	VCLatticeFigureOfMeritToCheckSymmetry(const Double& rhs);	// Sets only m_latfom.m_determ_GramMat = rhs;
	VCLatticeFigureOfMeritToCheckSymmetry(const BravaisType& ebrat,
										const SymMat43_VCData& S_red,
										const ePeakShiftFunctionType& type,
										const Double& wave_length,
										const vector<ZParawError>& peak_shift_param_rad);
	virtual ~VCLatticeFigureOfMeritToCheckSymmetry(){};
	
	void setLatticeFigureOfMerit(const LatticeFigureOfMeritZeroShift& arg) { m_latfom = arg; };
	const LatticeFigureOfMeritZeroShift& putLatticeFigureOfMerit() const { return m_latfom; };

	// Computes Mwu.
	inline void setWuFigureOfMerit(const Int4& num_ref_figure_of_merit,
									const vector<QData>& qdata,
									const Double& min_thred_num_hkl,
									const Double& max_thred_num_hkl) { m_latfom.setWuFigureOfMerit(num_ref_figure_of_merit, qdata, min_thred_num_hkl, max_thred_num_hkl); };

	// Sets the member m_latfom.m_figures_of_merit.
	inline void setFigureOfMerit(const Int4& num_ref_figure_of_merit,
									const vector<QData>& qdata,
									vector< VecDat3<Int4> >& closest_hkl_tray,
									Vec_BOOL& Q_observed_flag){ m_latfom.setFigureOfMerit(num_ref_figure_of_merit, qdata, closest_hkl_tray, Q_observed_flag); };

	// Returns true if the lattice has at least the symmetry of eps. 
	// On output, ans equals the equivalent lattice with symmetry of eps.
	bool checkLatticeSymmetry(const ePointGroup& epg_new, const Double& cv2,
								map< SymMat<VCData>, NRMat<Int4> >& ans) const;

	// Set-functions.
	inline void setLabel(const lattice_label& arg) { m_label = arg; };

	inline const lattice_label& putLatticeLabel() const { return m_label; };
	inline string putStringLabel() const { return CS_LABEL[(size_t)m_latfom.enumBravaisType()] + "0" + num2str<Int4>(m_label); };

	inline const SymMat43_VCData& putInitialForm() const { return m_S_red; };
	inline const SymMat<VCData> putInitialSellingReducedForm() const { return transform_sym_matrix(m_S_red.second, m_S_red.first); };

	inline Int4& putNumberOfLatticesInNeighborhood() { return m_num_lattice_found; };
	inline const Int4& putNumberOfLatticesInNeighborhood() const { return m_num_lattice_found; };

	void putLatticesOfHigherSymmetry(
			const ePointGroup& epg, const Double& cv2,
			vector<VCLatticeFigureOfMeritToCheckSymmetry>& lattice_result) const;

	void printLatticeInformation(const vector<VCLatticeFigureOfMeritToCheckSymmetry> lattice_result[],
						const OutputInfo outinfo[],
						const eABCaxis& abc_axis,
						const eRHaxis& rh_axis,
						const Double& resol,
						const Int4& label_start0,
						ostream* os) const;

	static bool cmpFOMdeWolff(const VCLatticeFigureOfMeritToCheckSymmetry& lhs, const VCLatticeFigureOfMeritToCheckSymmetry& rhs)
	{
		return LatticeFigureOfMerit::cmpFOMdeWolff( lhs.putLatticeFigureOfMerit(), rhs.putLatticeFigureOfMerit() );
	}
	static bool cmpFOMWu(const VCLatticeFigureOfMeritToCheckSymmetry& lhs, const VCLatticeFigureOfMeritToCheckSymmetry& rhs)
	{
		return LatticeFigureOfMerit::cmpFOMWu( lhs.putLatticeFigureOfMerit(), rhs.putLatticeFigureOfMerit() );
	}
	static bool cmpFOMReversed(const VCLatticeFigureOfMeritToCheckSymmetry& lhs, const VCLatticeFigureOfMeritToCheckSymmetry& rhs)
	{
		return LatticeFigureOfMerit::cmpFOMReversed( lhs.putLatticeFigureOfMerit(), rhs.putLatticeFigureOfMerit() );
	}
	static bool cmpFOMSymmetric(const VCLatticeFigureOfMeritToCheckSymmetry& lhs, const VCLatticeFigureOfMeritToCheckSymmetry& rhs)
	{
		return LatticeFigureOfMerit::cmpFOMSymmetric( lhs.putLatticeFigureOfMerit(), rhs.putLatticeFigureOfMerit() );
	}

	static bool cmpNumberOfNeighbors(const VCLatticeFigureOfMeritToCheckSymmetry& lhs, const VCLatticeFigureOfMeritToCheckSymmetry& rhs)
	{
		return lhs.putNumberOfLatticesInNeighborhood() > rhs.putNumberOfLatticesInNeighborhood();
	}

	// For GUI
	const lattice_label        &getref_m_label          () const {return m_label;}
	      lattice_label        &getref_m_label          ()       {return m_label;}
	const LatticeFigureOfMerit &getref_m_latfom         () const {return m_latfom;}
	      LatticeFigureOfMerit &getref_m_latfom         ()       {return m_latfom;}
	const SymMat43_VCData      &getref_m_S_red          () const {return m_S_red;}
	      SymMat43_VCData      &getref_m_S_red          ()       {return m_S_red;}
	const Int4                 &getref_num_lattice_found() const {return this->putNumberOfLatticesInNeighborhood();}
	      Int4                 &getref_num_lattice_found()       {return this->putNumberOfLatticesInNeighborhood();}
//	const vector<EquivLattice> *getref_lattice_equiv    () const {return m_lattice_equiv;}
//	      vector<EquivLattice> *getref_lattice_equiv    ()       {return m_lattice_equiv;}
};

inline bool operator<(const VCLatticeFigureOfMeritToCheckSymmetry& lhs, const VCLatticeFigureOfMeritToCheckSymmetry& rhs)
{
	return lhs.putLatticeFigureOfMerit() < rhs.putLatticeFigureOfMerit();
}

#endif /*VCLatticeFigureOfMeritToCheckSymmetry_HH_*/

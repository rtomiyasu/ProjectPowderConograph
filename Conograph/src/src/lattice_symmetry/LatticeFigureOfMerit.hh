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
#ifndef LATTICEFIGUREOFMERIT_HH_
#define LATTICEFIGUREOFMERIT_HH_

#include "../utility_data_structure/nrutil_nr.hh"
#include "../utility_data_structure/SymMatNplus1N.hh"
#include "../utility_data_structure/VecDat3.hh"
#include "../utility_func/transform_sym_matrix.hh"
#include "../utility_func/lattice_constant.hh"
#include "../utility_func/chToDouble.hh"
#include "../utility_lattice_reduction/matrix_NbyN.hh"
#include "../centring_type/enumCentringType.hh"
#include "../point_group/enumPointGroup.hh"
#include "../laue_group/LaueGroup.hh"
#include "../zparam/ZParawError.hh"
#include "../model_function/profile_function/global_function/enumPeakShiftFunctionType.hh"
#include "../bravais_type/BravaisType.hh"
#include "enumSortCriterion.hh"
#include "HKL_Q.hh"
#include "check_equiv.hh"


// Class for outputting information about a lattice in index file.
class LatticeFigureOfMerit
{
	friend inline bool operator<(const LatticeFigureOfMerit& lhs, const LatticeFigureOfMerit& rhs);

public:
	class SetOfFigureOfMerit
	{
	private:
		Int4 m_num_ref_figure_of_merit;
		Int4 m_num_q_observed;
		Double m_num_total_hkl;
		Double m_figure_of_merit_Wolff;
		Double m_figure_of_merit_Wu;
		Double m_reversed_figure_of_merit;

		static string putStrFigureOfMeritWolff(const Int4& num_ref_fom){ return putLabel(SCM) + num2str(num_ref_fom); };
		static string putStrFigureOfMeritWu(const Int4& num_ref_fom){ return putLabel(SCWuM) + num2str(num_ref_fom); };
		static string putStrReversedFigureOfMeritWolff(const Int4& num_ref_fom){ return putLabel(SCRevM) + num2str(num_ref_fom); };
		static string putStrSymmetricFigureOfMeritWolff(const Int4& num_ref_fom){ return putLabel(SCSymM) + num2str(num_ref_fom); };

	public:
		SetOfFigureOfMerit(){ this->reset(); }
		~SetOfFigureOfMerit(){}

		inline void reset()
		{
			m_num_ref_figure_of_merit = 0;
			m_num_q_observed = 0;
			m_num_total_hkl = 0.0;
			m_figure_of_merit_Wolff = 0.0;
			m_figure_of_merit_Wu = 0.0;
			m_reversed_figure_of_merit = 0.0;
		}

		inline Int4& putNumberOfReflectionsForFigureOfMerit() { return m_num_ref_figure_of_merit; };
		inline Int4& putNumQobsAssociatedWithCloseHKL() { return m_num_q_observed; };
		inline Double& putContinuousNumberOfHKLInRange() { return m_num_total_hkl; };

		inline const Int4& putNumQobsAssociatedWithCloseHKL() const { return m_num_q_observed; };
		inline const Int4& putNumberOfReflectionsForFigureOfMerit() const { return m_num_ref_figure_of_merit; };
		inline const Double& putContinuousNumberOfHKLInRange() const { return m_num_total_hkl; };

		inline Double& putFigureOfMeritWolff() { return m_figure_of_merit_Wolff; };
		inline Double& putFigureOfMeritWu() { return m_figure_of_merit_Wu; };
		inline Double& putReversedFigureOfMerit() { return m_reversed_figure_of_merit; };

		inline const Double& putFigureOfMeritWolff() const { return m_figure_of_merit_Wolff; };
		inline const Double& putFigureOfMeritWu() const { return m_figure_of_merit_Wu; };
		inline const Double& putReversedFigureOfMerit() const { return m_reversed_figure_of_merit; };

		inline string putLabel_FigureOfMeritWolff() const { return putStrFigureOfMeritWolff(m_num_ref_figure_of_merit); };
		inline string putLabel_FigureOfMeritWu() const { return putStrFigureOfMeritWu(m_num_ref_figure_of_merit); };
		inline string putLabel_ReversedFigureOfMeritWolff() const { return putStrReversedFigureOfMeritWolff(m_num_ref_figure_of_merit); };
		inline string putLabel_SymmetricFigureOfMeritWolff() const { return putStrSymmetricFigureOfMeritWolff(m_num_ref_figure_of_merit); };

		inline Double putSymmetricFigureOfMerit() const { return m_reversed_figure_of_merit * m_figure_of_merit_Wolff; };
	};

private:
	static const NRMat<Int4> m_tmat_prim_to_face;
	static const NRMat<Int4> m_tmat_prim_to_body;
	static const NRMat<Int4> m_tmat_prim_to_rhomhex;
	static const NRMat<Int4> m_tmat_prim_to_base[3];
	static const NRMat<Int4> m_tmat_prim_to_prim;

	BravaisType m_brat;

	// This matrix is a result of optimization, therefore, it may be not reduced.
	// m_S_optimized.second * m_S_optimized.first * Transpose(m_S_optimized.second) is obtuse.
	SymMat43_Double m_S_optimized;

	// The inverse of m_S_red is Buerger-reduced.
	SymMat<Double> m_S_red;

	Double m_determ_S_red;
	
	SetOfFigureOfMerit m_figures_of_merit;
	
	// Sets m_S_red from m_S_optimized.
	// On output, trans_mat gives the matrix such that trans_mat * m_S_red * transpose(trans_mat) = m_S_optimized.first.
	void setInverseOfBuergerReducedForm(NRMat<Int4>& trans_mat);

	// Returns unit-cell parameters with other centrings using m_S_red.
	void putEquivalentLatticeConstantsDegreeWithOtherCentring(const eABCaxis& abc_axis,
																const eRHaxis& rh_axis,
																const Double& resol,
																vector< pair< eBravaisType, SymMat<Double> > >& ans) const;

protected:
	static const Double m_cv2;

public:
	LatticeFigureOfMerit();
	LatticeFigureOfMerit(const Double& rhs);	// Sets only m_determ_GramMat = rhs;
	LatticeFigureOfMerit(const BravaisType& ebrat,
							const SymMat43_Double& S_red);
	virtual ~LatticeFigureOfMerit(){}; 
	
	void putMillerIndicesInRange(const Double& qrange_end,
									const Int4& irc_type,
									vector<HKL_Q>& cal_hkl_tray) const;

	void setDeWolffFigureOfMerit(const Int4& num_ref_figure_of_merit, const Int4& irc_type, const vector<QData>& qdata);

	void setWuFigureOfMerit(const Int4& num_ref_figure_of_merit, const vector<QData>& qdata,
					const Double& min_thred_num_hkl,
					const Double& max_thred_num_hkl);

	// Return false if Qdata is not set or
	// the number of unindexed reflections is larger max_num_false_peak.
	void setFigureOfMerit(const Int4& num_ref_figure_of_merit, const vector<QData>& qdata,
							vector< VecDat3<Int4> >& closest_hkl_tray,
							Vec_BOOL& Q_observed_flag);

	inline void setFigureOfMerit(const Int4& num_ref_figure_of_merit, const vector<QData>& qdata)
	{
		vector< VecDat3<Int4> > closest_hkl_tray;
		Vec_BOOL Q_observed_flag;
		setFigureOfMerit(num_ref_figure_of_merit, qdata, closest_hkl_tray, Q_observed_flag);
	};

	//  true : NormM has been improved.
	// false : NormM has not been improved.
	pair<bool, ZErrorMessage> fitLatticeParameterLinear(const vector<QData>& qdata,
													const vector< VecDat3<Int4> >& hkl_to_fit,
													const vector<bool>& fix_or_fit_flag,
													const bool& output_view_flag);

	// Change the lattice constants to string.
	inline string printOptimizedLatticeConstants(const eABCaxis& axis1,
													const eRHaxis& axis2,
													const Int4& precision) const;

	// Output information on the lattice.
	void printLatticeInformation(const eABCaxis& abc_axis,
						const eRHaxis& rh_axis,
						const Double& resol,
						const Int4& label_start0,
						ostream* os) const;

	// Set-functions.
	// This method assumes that S.second * S.first * Transpose(S.second) is obtuse.
	void setLatticeConstants43(const BravaisType& brat, const SymMat43_Double& S);

	ZErrorMessage setLatticeConstants(const BravaisType& brat, const SymMat<Double>& S);

	Int4 checkDominantZone() const;

	// Replace m_S_optimized by m_S_red.
	// On output, trans_mat gives the matrix such that trans_mat * m_S_red * transpose(trans_mat) equals the original m_S_optimized.
	inline void reduceLatticeConstants(NRMat<Int4>& trans_mat);

	// The unit of alpha, beta, gamma is degree.
	inline ZErrorMessage setLatticeConstantsDegree(const BravaisType& brat,
													const VecDat3<Double>& length_axis,
													const VecDat3<Double>& angle_axis);

	static void putLatticeConstantsDegree(const BravaisType& brat,
													const SymMat<Double>& S,
													const eABCaxis& axis1,
													const eRHaxis& axis2,
													VecDat3<Double>& length, VecDat3<Double>& angle);

	// Put-functions (Returns a value.)
	inline void putOptimizedLatticeConstantsDegree(const eABCaxis& axis1,
													const eRHaxis& axis2,
													VecDat3<Double>& length, VecDat3<Double>& angle) const;

	inline void putReducedLatticeConstantsDegree(const eABCaxis& axis1,
									const eRHaxis& axis2,
									VecDat3<Double>& length, VecDat3<Double>& angle) const;

	inline Double putLatticeVolume() const { return 1.0 / sqrt(m_determ_S_red); };

	// Put-functions (Returns a constant reference.)
	inline const Double& putCriticalValueSquare() const { return m_cv2; };

	inline const ePointGroup& enumLaueGroup() const { return m_brat.enumLaueGroup(); };
	inline const eCentringType& enumCentringType() const { return m_brat.enumCentringType(); };
	inline const eBravaisType& enumBravaisType() const { return m_brat.enumBravaisType(); };
	inline const BravaisType& putBravaisType() const { return m_brat; };

	// The returned matrix is Buerger-reduced matrix equivalent with m_S_red_optimized.
	inline const SymMat<Double>& putInverseOfBuergerReducedForm() const { return m_S_red; };

	inline const SymMat43_Double& putOptimizedForm() const { return m_S_optimized; };
	inline SymMat<Double> putSellingReducedForm() const { return transform_sym_matrix(m_S_optimized.second, m_S_optimized.first); };

	inline const Double& putDeterminantOfGramMatrix() const { return m_determ_S_red; };
	inline const SetOfFigureOfMerit& putFiguresOfMerit() const { return m_figures_of_merit; };

	// Put-functions (Returns a non-constant reference.)
	// Returns almost equivalent unit-cell parameters of different centring-types.
	void putEquivalentLatticeConstantsDegreeWithOtherCentring(const eABCaxis& abc_axis,
																const eRHaxis& rh_axis,
																const Double& resol,
																vector< pair< eBravaisType, pair< VecDat3<Double>, VecDat3<Double> > > >& ans) const;

	static bool cmpFOMdeWolff(const LatticeFigureOfMerit& lhs, const LatticeFigureOfMerit& rhs)
	{
		return lhs.m_figures_of_merit.putFigureOfMeritWolff() > rhs.m_figures_of_merit.putFigureOfMeritWolff();
	}

	static bool cmpFOMWu(const LatticeFigureOfMerit& lhs, const LatticeFigureOfMerit& rhs)
	{
		return lhs.m_figures_of_merit.putFigureOfMeritWu() > rhs.m_figures_of_merit.putFigureOfMeritWu();
	}

	static bool cmpFOMReversed(const LatticeFigureOfMerit& lhs, const LatticeFigureOfMerit& rhs)
	{
		return lhs.m_figures_of_merit.putReversedFigureOfMerit() > rhs.m_figures_of_merit.putReversedFigureOfMerit();
	}

	static bool cmpFOMSymmetric(const LatticeFigureOfMerit& lhs, const LatticeFigureOfMerit& rhs)
	{
		return lhs.m_figures_of_merit.putSymmetricFigureOfMerit() > rhs.m_figures_of_merit.putSymmetricFigureOfMerit();
	}

	// For GUI
	const BravaisType            &getref_m_brat()             const {return m_brat;}
	      BravaisType            &getref_m_brat()                   {return m_brat;} 
	const SetOfFigureOfMerit     &getref_m_figures_of_merit() const {return m_figures_of_merit;}
	      SetOfFigureOfMerit     &getref_m_figures_of_merit()       {return m_figures_of_merit;}
	const SymMat43_Double        &getref_m_S_red_optimized()  const {return m_S_optimized;}
	      SymMat43_Double        &getref_m_S_red_optimized()        {return m_S_optimized;}
	const Double                 &getref_m_determ_S_red()     const {return m_determ_S_red;}
	      Double                 &getref_m_determ_S_red()           {return m_determ_S_red;}
	const SymMat<Double>     &getref_m_S_red()            const {return m_S_red;}
	      SymMat<Double>     &getref_m_S_red()                  {return m_S_red;}

};


inline string LatticeFigureOfMerit::printOptimizedLatticeConstants(const eABCaxis& axis1,
		const eRHaxis& axis2,
		const Int4& precision) const
{
	VecDat3<Double> length_axis, angle_axis;
	putOptimizedLatticeConstantsDegree(axis1, axis2, length_axis, angle_axis);
	return chToString(length_axis, angle_axis, precision);
}


inline void LatticeFigureOfMerit::putOptimizedLatticeConstantsDegree(const eABCaxis& axis1,
		const eRHaxis& axis2,
		VecDat3<Double>& length_axis, VecDat3<Double>& angle_axis) const
{
	putLatticeConstantsDegree( m_brat, m_S_optimized.first, axis1, axis2, length_axis, angle_axis );
}


inline void LatticeFigureOfMerit::putReducedLatticeConstantsDegree(
		const eABCaxis& axis1,
		const eRHaxis& axis2, VecDat3<Double>& length_axis, VecDat3<Double>& angle_axis) const
{
	putLatticeConstantsDegree( m_brat, this->putInverseOfBuergerReducedForm(),
								axis1, axis2, length_axis, angle_axis );
}



inline ZErrorMessage LatticeFigureOfMerit::setLatticeConstantsDegree(const BravaisType& brat,
		const VecDat3<Double>& length,
		const VecDat3<Double>& angle)
{
	LaueGroup lg(brat.enumLaueGroup());
	ZErrorMessage zerr = lg->checkLatticeConstantError(length, angle);
	if( zerr.putErrorType() != ZErrorNoError ) return zerr;

	SymMat<Double> Sval(3);
	calCoParameter(length, angle, Sval);

	return this->setLatticeConstants(brat, Sval);
}

inline void LatticeFigureOfMerit::reduceLatticeConstants(NRMat<Int4>& trans_mat)
{
	// trans_mat * m_S_red * transpose(trans_mat) = m_S_optimized.first(old).
	setInverseOfBuergerReducedForm(trans_mat);
	// m_S_optimized.first(new) = m_S_red.
	m_S_optimized.first = m_S_red;

	// m_S_optimized.second(new) * m_S_optimized.first(new) * transpose(m_S_optimized.second(new))
	//    = m_S_optimized.second(old) * m_S_optimized.first(old) * transpose(m_S_optimized.second(old)) is Delone reduced.
	m_S_optimized.second = mprod(m_S_optimized.second, trans_mat);
}


template<class T>
inline T norm(const VecDat3<Int4>& lhs, const SymMat<T>& rhs)
{
	assert( rhs.size() >= 3 );
	
	return rhs(0,0)*(lhs[0]*lhs[0]) + rhs(1,1)*(lhs[1]*lhs[1]) + rhs(2,2)*(lhs[2]*lhs[2])
			+ ( rhs(0,1)*(lhs[0]*lhs[1]) + rhs(0,2)*(lhs[0]*lhs[2]) + rhs(1,2)*(lhs[1]*lhs[2]) )*2.0;
}


inline bool operator<(const LatticeFigureOfMerit& lhs, const LatticeFigureOfMerit& rhs)
{ 
	return lhs.m_determ_S_red < rhs.m_determ_S_red;
}

void putTransformMatrixToBuergerReduced(const SymMat<Double>& S, NRMat<Int4>& trans_mat);

#endif /*LATTICEFIGUREOFMERIT_HH_*/

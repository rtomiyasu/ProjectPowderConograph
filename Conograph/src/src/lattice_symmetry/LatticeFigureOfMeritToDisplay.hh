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
#ifndef LatticeFigureOfMeritToDisplay_HH_
#define LatticeFigureOfMeritToDisplay_HH_

#include <map>
#include "LatticeFigureOfMeritZeroShift.hh"
#include "../utility_data_structure/range.hh"
#include "../qc/reflection_conditions.hh"


// Class for outputting information about a lattice in IGOR file.
class LatticeFigureOfMeritToDisplay
{
private:
	LatticeFigureOfMeritZeroShift m_latfom;
	
	vector<QData> m_qdata;
	vector< multimap<Double, vector<HKL_Q>::const_iterator> > m_associated_hkl_tray;
	vector<HKL_Q> m_cal_hkl_tray;
	
	// The size of the following arrays equals zero or m_figures_of_merit.m_num_ref_figure_of_merit.
	vector< VecDat3<Int4> > m_hkl_to_fit;
	vector<bool> m_fix_or_fit_flag;	// 0:fix, 1:fit.

	// If this is negative value, only the reflection conditions derived from the Bravais type are considered.
	Int4 m_type_of_reflection_conditions;

	// for GUI
	bool m_showsTicks;

public:
	void putPeakPosToFit(const ControlParam& cdata, Vec_DP& cal_q_tray, Vec_DP& cal_pos_tray) const;
	void putStandardMillerIndicesToFit(vector< VecDat3<Int4> >& hkl_to_fit,
										vector<bool>& fix_or_fit_flag) const;

	LatticeFigureOfMeritToDisplay();
	~LatticeFigureOfMeritToDisplay(){};

	inline const eBravaisType& enumBravaisType() const { return m_latfom.enumBravaisType(); };
	inline const LatticeFigureOfMeritZeroShift& putLatticeFigureOfMerit() const { return m_latfom; };
	inline const vector<HKL_Q>& putCalMillerIndices() const { return m_cal_hkl_tray; };
	inline const vector< multimap<Double, vector<HKL_Q>::const_iterator> >& putAssociatedMillerIndices() const { return m_associated_hkl_tray; };

	// axis1: axis of the argument lattice constants for Monoclinic(C) or Orthorhombic(C).
	// axis2: axis of the argument lattice constants for Rhombohedral.
	inline ZErrorMessage setLatticeConstantsDegree(const eBravaisType& ecys,
													const eABCaxis& axis1,
													const eRHaxis& axis2,
													const VecDat3<Double>& length_axis,
													const VecDat3<Double>& angle_axis) { return m_latfom.setLatticeConstantsDegree(BravaisType(ecys, axis1, axis2), length_axis, angle_axis); };
	inline ZErrorMessage setPeakShiftParamDegree(const ePeakShiftFunctionType& type,
													const Double& wave_length,
													const vector<ZParawError>& peak_shift_param_deg,	// The errors are not used in this method.
													const PeakPosData& pdata);
	inline void reduceLatticeConstants();

	inline void putOptimizedLatticeConstantsDegree(const eABCaxis& axis1,
													const eRHaxis& axis2,
													VecDat3<Double>& length_axis,
													VecDat3<Double>& angle_axis) { m_latfom.putOptimizedLatticeConstantsDegree(axis1, axis2, length_axis, angle_axis); };

	// axis1: axis of the argument lattice constants for Monoclinic(C) or Orthorhombic(C).
	// axis2: axis of the argument lattice constants for Rhombohedral.
	inline void putReducedLatticeConstantsDegree(const eABCaxis& axis1,
													const eRHaxis& axis2,
													VecDat3<Double>& length_axis,
													VecDat3<Double>& angle_axis) { m_latfom.putReducedLatticeConstantsDegree(axis1, axis2, length_axis, angle_axis); };

	inline const Double& putWaveLength() const { return m_latfom.putWaveLength(); };
	inline const ePeakShiftFunctionType& putPeakShiftFunctionType() const { return m_latfom.putPeakShiftFunctionType(); };
	inline vector<ZParawError> putPeakShiftParamDegree() const { return m_latfom.putPeakShiftParamDegree(); };

	inline const vector<QData>& putQDataModifiedWithNewPeakShiftParam() const { return m_qdata; };

	inline void setDeWolffFigureOfMerit(const Int4& num_ref_figure_of_merit){ m_latfom.setDeWolffFigureOfMerit(num_ref_figure_of_merit, m_type_of_reflection_conditions, m_qdata); };
	inline void setFigureOfMerit(const Int4& num_ref_figure_of_merit){ m_latfom.setFigureOfMerit(num_ref_figure_of_merit, m_qdata); };

	void setTypeOfSystematicAbsences(const Int4& arg);
	inline string putShortStringTypeOfSystematicAbsences() const
	{
		return putInformationOnReflectionConditions(this->putLatticeFigureOfMerit().putBravaisType(), m_type_of_reflection_conditions).putShortStringType();
	};
	inline const string& putStringTypeOfSystematicAbsences() const
	{
		return putInformationOnReflectionConditions(this->putLatticeFigureOfMerit().putBravaisType(), m_type_of_reflection_conditions).putStringType();
	};
	inline const string& putStringReflectionConditions() const
	{
	     return putInformationOnReflectionConditions(this->putLatticeFigureOfMerit().putBravaisType(), m_type_of_reflection_conditions).putStringConditions();
	};
	inline const DataReflectionConditions& putDataOnReflectionConditions() const
	{
	     return putInformationOnReflectionConditions(this->putLatticeFigureOfMerit().putBravaisType(), m_type_of_reflection_conditions);
	};


	// Resets m_associated_hkl_tray and q-values in m_cal_hkl_tray.
	void resetQValuesInRange(const Int4& num_fit_data);
	// Resets m_associated_hkl_tray and m_cal_hkl_tray.
	void resetMillerIndicesInRange(const Int4& num_fit_data);
	// Resets m_hkl_to_fit and m_fix_or_fit_flag.
	void resetMillerIndicesToFit() { putStandardMillerIndicesToFit(m_hkl_to_fit, m_fix_or_fit_flag); };

	inline ZErrorMessage setFittingIDs(const vector<bool>& arg);
	inline const vector<bool>& putFittingIDs() const { return m_fix_or_fit_flag; };
	inline ZErrorMessage setMillerIndicesToFit(const vector< VecDat3<Int4> >& arg);
	inline const vector< VecDat3<Int4> >& putMillerIndicesToFit() const { return m_hkl_to_fit; };

	// Return peak positions for this lattice candidate.
	void putCalculatedPeakPosInRange(const ControlParam& cdata, Vec_DP& cal_pos_tray) const;

	// Returned value
	// > 1 : Optimization has succeeded at least twice. (The members m_hkl_to_fit and m_fix_or_fit_flag are changed.)
	// > 0 : Optimization has succeeded. (Lattice constants are changed.)
	// 0 : Optimization has failed.
	Int4 fitLatticeParameter(const PeakPosData& pdata, const vector<etype_ID>& fitflag,
								const Int4& Max_ITNUM,
								const Double& limiter);

	// Output indexing results.
	void printIndexingResult(const ControlParam& cdata,
						const PeakPosData& pdata,
						const Int4& label_start0,
						ostream* os) const;

	// For GUI.
	const vector<bool>           &getref_m_fix_or_fit_flag()  const {return m_fix_or_fit_flag;}
	      vector<bool>           &getref_m_fix_or_fit_flag()        {return m_fix_or_fit_flag;}
	const vector<VecDat3<Int4> > &getref_m_hkl_to_fit()       const {return m_hkl_to_fit;}
	      vector<VecDat3<Int4> > &getref_m_hkl_to_fit()             {return m_hkl_to_fit;}
	const LatticeFigureOfMeritZeroShift                             &getref_m_latfom()              const {return m_latfom;}
	      LatticeFigureOfMeritZeroShift                             &getref_m_latfom()                    {return m_latfom;}
	const vector<QData>                                             &getref_m_qdata()               const {return m_qdata;}
	      vector<QData>                                             &getref_m_qdata()                     {return m_qdata;}
	const vector< multimap<Double, vector<HKL_Q>::const_iterator> > &getref_m_associated_hkl_tray() const {return m_associated_hkl_tray;}
	      vector< multimap<Double, vector<HKL_Q>::const_iterator> > &getref_m_associated_hkl_tray()       {return m_associated_hkl_tray;}
	const vector<HKL_Q>                                             &getref_m_cal_hkl_tray()        const {return m_cal_hkl_tray;}
	      vector<HKL_Q>                                             &getref_m_cal_hkl_tray()              {return m_cal_hkl_tray;}
  	const bool                                                      &getref_m_showsTicks()          const {return m_showsTicks;}
  	      bool                                                      &getref_m_showsTicks()                {return m_showsTicks;}
};


inline ZErrorMessage LatticeFigureOfMeritToDisplay::setFittingIDs(const vector<bool>& arg)
{
	if( arg.size() != m_associated_hkl_tray.size() )
	{
		return ZErrorMessage(ZErrorArgmentSize, __FILE__, __LINE__, __FUNCTION__);
	}
	 m_fix_or_fit_flag = arg;
	return ZErrorMessage();
}

inline ZErrorMessage LatticeFigureOfMeritToDisplay::setMillerIndicesToFit(const vector< VecDat3<Int4> >& arg)
{
	if( arg.size() != m_associated_hkl_tray.size() )
	{
		return ZErrorMessage(ZErrorArgmentSize, __FILE__, __LINE__, __FUNCTION__);
	}
	m_hkl_to_fit = arg;
	return ZErrorMessage();
};


inline ZErrorMessage LatticeFigureOfMeritToDisplay::setPeakShiftParamDegree(
		const ePeakShiftFunctionType& type,
		const Double& wave_length,
		const vector<ZParawError>& peak_shift_param_deg,
		const PeakPosData& pdata)
{
	static const Double RadDeg =  PI() / 180.0;

	vector<ZParawError> peak_shift_param_rad = peak_shift_param_deg;
	for(vector<ZParawError>::iterator it=peak_shift_param_rad.begin(); it<peak_shift_param_rad.end(); it++)
	{
		*it *= RadDeg;
	}
	return m_latfom.setPeakShiftParamRadian(VCData::putPeakQData(), type, wave_length, peak_shift_param_rad, pdata,
											VCData::putPeakQData().size(), m_qdata);
};


inline VecDat3<Int4> product_hkl(const VecDat3<Int4>& lhs, const NRMat<Int4>& rhs)
{
	assert( rhs.nrows() >= 3 && rhs.ncols() == 3 );

	VecDat3<Int4> ans;
	ans[0] = lhs[0]*rhs[0][0] + lhs[1]*rhs[1][0] + lhs[2]*rhs[2][0];
	ans[1] = lhs[0]*rhs[0][1] + lhs[1]*rhs[1][1] + lhs[2]*rhs[2][1];
	ans[2] = lhs[0]*rhs[0][2] + lhs[1]*rhs[1][2] + lhs[2]*rhs[2][2];
	return ans;
}


inline void LatticeFigureOfMeritToDisplay::reduceLatticeConstants()
{
	NRMat<Int4> trans_mat;
	// trans_mat * m_latfom.m_S_optimized.first(new) * transpose(trans_mat) = m_latfom.m_S_optimized.first(old).
	m_latfom.reduceLatticeConstants(trans_mat);

	for(vector<HKL_Q>::iterator it=m_cal_hkl_tray.begin(); it<m_cal_hkl_tray.end(); it++)
	{
		it->setHKL( product_hkl(it->HKL(), trans_mat) );
	}

	for(vector< VecDat3<Int4> >::iterator it=m_hkl_to_fit.begin(); it<m_hkl_to_fit.end(); it++)
	{
		*it = product_hkl(*it, trans_mat);
	}
};

inline bool cmpDeWolff(const LatticeFigureOfMeritToDisplay& lhs, const LatticeFigureOfMeritToDisplay& rhs)
{
	return lhs.putLatticeFigureOfMerit().putFiguresOfMerit().putFigureOfMeritWolff()
			> rhs.putLatticeFigureOfMerit().putFiguresOfMerit().putFigureOfMeritWolff();
}


#endif /*LatticeFigureOfMeritToDisplay_HH_*/

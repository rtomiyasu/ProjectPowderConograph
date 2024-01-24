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
#ifndef LATTICEFIGUREOFMERITZEROSHIFT_HH_
#define LATTICEFIGUREOFMERITZEROSHIFT_HH_

#include "LatticeFigureOfMerit.hh"

class ControlParam;
class PeakPosData;

// Class for outputting information about a lattice in index file.
class LatticeFigureOfMeritZeroShift : public LatticeFigureOfMerit
{
private:
	ePeakShiftFunctionType m_etype_peak_shift;
	Double m_WlengthX;
	vector<ZParawError> m_peak_shift_param_rad;

public:
	LatticeFigureOfMeritZeroShift();
	LatticeFigureOfMeritZeroShift(const Double& rhs);	// Sets only m_determ_GramMat = rhs;
	virtual ~LatticeFigureOfMeritZeroShift(){};
	
	// On output, the size of qdata_modified equals qdata_modified_size.
	//  true : NormM has been improved.
	//  false : NormM has not been improved.
	pair<bool, ZErrorMessage> fitLatticeParameter(
							const vector<QData>& qdata,
							const vector< VecDat3<Int4> >& hkl_to_fit,
							const vector<bool>& fix_or_fit_flag,
							const bool& output_view_flag,
							const PeakPosData& pdata,
							const vector<etype_ID>& peak_shift_fitflag,
							const Int4& Max_ITNUM,
							const Int4& qdata_modified_size,
							vector<QData>& qdata_modified);

	inline const ePeakShiftFunctionType& putPeakShiftFunctionType() const { return m_etype_peak_shift; };
	inline const Double& putWaveLength() const { return m_WlengthX; };
	inline const vector<ZParawError>& putPeakShiftParamRadian() const { return m_peak_shift_param_rad; };
	inline vector<ZParawError> putPeakShiftParamDegree() const;

	// On output, the size of qdata_modified equals qdata_modified_size.
	ZErrorMessage setPeakShiftParamRadian(const ePeakShiftFunctionType& type,
							const Double& wave_length,
							const vector<ZParawError>& peak_shift_param_rad);

	// On output, the size of qdata_modified equals qdata_modified_size.
	ZErrorMessage setPeakShiftParamRadian(const vector<QData>& qdata,
							const ePeakShiftFunctionType& type,
							const Double& wave_length,
							const vector<ZParawError>& peak_shift_param_error,
							const PeakPosData& pdata,
							const Int4& qdata_modified_size,
							vector<QData>& qdata_modified);

	// For GUI.
	const ePeakShiftFunctionType &getref_m_etype_peak_shift()      const {return m_etype_peak_shift;}
	      ePeakShiftFunctionType &getref_m_etype_peak_shift()            {return m_etype_peak_shift;}
	const Double                 &getref_m_WlengthX()              const {return m_WlengthX;}
	      Double                 &getref_m_WlengthX()                    {return m_WlengthX;}
	const vector<ZParawError>    &getref_m_peak_shift_param_rad()  const {return m_peak_shift_param_rad;}
	      vector<ZParawError>    &getref_m_peak_shift_param_rad()        {return m_peak_shift_param_rad;}
	const LatticeFigureOfMerit   &getref_base()                    const {return *this;}
	      LatticeFigureOfMerit   &getref_base()                          {return *this;}
};

inline vector<ZParawError> LatticeFigureOfMeritZeroShift::putPeakShiftParamDegree() const
{
	static const Double DegRad =  180.0 / PI();

	vector<ZParawError> peak_shift_deg = m_peak_shift_param_rad;
	for(vector<ZParawError>::iterator it=peak_shift_deg.begin(); it<peak_shift_deg.end(); it++)
	{
		*it *= DegRad;
	}
	return peak_shift_deg;
}

#endif /*LATTICEFIGUREOFMERITZEROSHIFT_HH_*/

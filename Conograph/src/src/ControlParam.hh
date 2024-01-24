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
#ifndef _ControlParam_h_
#define _ControlParam_h_
// ControlParam.hh

#include "RietveldAnalysisTypes.hh"
#include "zerror_type/error_out.hh"
#include "model_function/profile_function/global_function/enumPeakShiftFunctionType.hh"
#include "bravais_type/enumBravaisType.hh"
#include "utility_func/zmath.hh"
#include "utility_rw_param/I_ReadData.hh"

class QData;

class ControlParam : public I_ReadData
{
private:
	enum{ NUM_LS = 14 };
	static const map<eABCaxis, string> ABCaxisLabel;
	static const map<eRHaxis, string> RHaxisLabel;
	
	// Parameters for search.
	static const pair< RWParamProperty, RWParamData<bool> > IsAngleDispersion_Data;
	static const pair< RWParamProperty, RWParamData<bool> > OutputSymmetry_Data[NUM_LS];
	static const pair< RWParamProperty, RWParamData<Vec_DP> > ConvParam_Data;
	static const pair< RWParamProperty, RWParamData<Vec_DP> > PeakShiftParam_Data;
	static const pair< RWParamProperty, RWParamData<Double> > WlengthX_Data;

	static const pair< RWParamProperty, RWParamData<Double> > CVForLinearSum_Data;	// critical value;

	static const pair< RWParamProperty, RWParamData<string> > str_MaxPeakNum_Data;
	static const pair< RWParamProperty, RWParamData<string> > str_MaxEdgeNum_Data;
	static const pair< RWParamProperty, RWParamData<string> > str_MaxNodeNum_Data;
	static const pair< RWParamProperty, RWParamData<string> > str_MinUnitCellVolume_Data;
	static const pair< RWParamProperty, RWParamData<string> > str_MaxUnitCellVolume_Data;

	// Parameters for output.
	static const pair< RWParamProperty, RWParamData<Int4> > NumRefFigureOfMerit_Data;

	static const pair< RWParamProperty, RWParamData<Double> > MinUnitCellEdgeABC_Data;
	static const pair< RWParamProperty, RWParamData<Double> > MaxUnitCellEdgeABC_Data;
	static const pair< RWParamProperty, RWParamData<Double> > MinFOM_Data;
	static const pair< RWParamProperty, RWParamData<Double> > Resol_Data;	// unit->Angstrom.
		
	static const pair< RWParamProperty, RWParamData<string> > str_MaxNumPeakInRange_Data;
	static const pair< RWParamProperty, RWParamData<string> > str_MinNumPeakInRange_Data;

	static const pair< RWParamProperty, RWParamData<Int4> > SearchLevel_Data;
	static const pair< RWParamProperty, RWParamData<Int4> > MaxSizeForEachBRAVAIS_Data;
	static const pair< RWParamProperty, RWParamData<string> > MonoBaseAxis_Data;
	static const pair< RWParamProperty, RWParamData<string> > RhomAxis_Data;
	static const pair< RWParamProperty, RWParamData<Double> > ThresholdNormM_Data;
	static const pair< RWParamProperty, RWParamData<Double> > ThresholdRevM_Data;
	static const pair< RWParamProperty, RWParamData<Double> > MinLatticePointDistance_Data;
	static const pair< RWParamProperty, RWParamData<Double> > ResoltoCheckComputedLines_Data;
	static const pair< RWParamProperty, RWParamData<Double> > MaxQtoCheckComputedLines_Data;
	const pair< RWParamProperty, RWParamData<Int4> > NumCores_Data;

	bool ReadConfigurationParameters;

	bool IsADorTOF;	// (0:tof, 1;angle dispersion)
	bool OutputSymmetry[NUM_LS];

	vector<Double> ConvParam;
	vector<Double> PeakShiftParam;	// Peak-shift parameters(rad.)
	Double WlengthX;    // If TOF, empty. If angle dispersion, the sum is always 1.
    
	// Parameters for search.
	Int4 SearchLevel;
	Int4 MaxPeakNum;
	Int4 MaxEdgeNum;
	Int4 MaxNodeNum;

	Double MinUnitCellVolume;
	Double MaxUnitCellVolume;
	Double CVForLinearSum;	// critical value;

	// Parameters for output.
	Int4 NumRefFigureOfMerit;
	Int4 MaxNumPeakInRange;
	Int4 MinNumPeakInRange;

	Double MinUnitCellEdgeABC;
	Double MaxUnitCellEdgeABC;
	Double MinFOM;
	Double Resol;	// unit->Angstrom.

	// Enviromental parameters.
	Int4 NumCores;
	Int4 MaxSizeForEachBRAVAIS;
	string MonoBaseAxis;
	string RhomAxis;
	Double ThresholdNormM;
	Double ThresholdRevM;
	Double MinLatticePointDistance;
	Double ResoltoCheckComputedLines;
	Double MaxQtoCheckComputedLines;	// unit: angstrom^{-2}.

	string str_MaxPeakNum;
	string str_MaxEdgeNum;
	string str_MaxNodeNum;
	string str_MinUnitCellVolume;
	string str_MaxUnitCellVolume;
	string str_MaxNumPeakInRange;
	string str_MinNumPeakInRange;

protected:
	ZErrorMessage checkData(const RWParam_void& param) const;
	ZErrorMessage checkIfDataAreSet(const RWParam_void& param, const Int4& num) const;

public:
    ControlParam();
    virtual ~ControlParam();

    inline void setReadConfigurationParameters(const bool& arg) { ReadConfigurationParameters = arg; };
    inline const bool& putReadConfigurationParameters() const { return ReadConfigurationParameters; };

    // Set functions.
    // Parameters for search.
    inline void setIsAngleDispersion(const bool& arg){ IsADorTOF = arg; };
	inline void setOutputSymmetry(const bool* arg){ for(Int4 i=0; i<NUM_LS; i++) OutputSymmetry[i] = arg[i]; };

	inline void putConvParam(const Vec_DP& arg){ ConvParam = arg; };
	inline void putPeakShiftParam(const Vec_DP& arg){ PeakShiftParam = arg; };
	inline void putWaveLength(const Double& arg){ WlengthX = arg; };

	inline void setMaxPeakNum(const Int4& arg) { MaxPeakNum = arg; };
	inline void setMaxEdgeNum(const Int4& arg) { MaxEdgeNum = arg; };
	inline void setMaxNodeNum(const Int4& arg) { MaxNodeNum = arg; };
	inline void setMaxUnitCellVolume(const Double& arg) { MaxUnitCellVolume = arg; };
	inline void setMinUnitCellVolume(const Double& arg) { MinUnitCellVolume = arg; };
	inline void setCriticalValueForLinearSum(const Double& arg) { CVForLinearSum = arg; };

    // Parameters for output.
	inline void setNumberOfReflectionsForFigureOfMerit(const Int4& arg) { NumRefFigureOfMerit = arg; };
	inline void setMinUnitCellEdgeABC(const Double& arg) { MinUnitCellEdgeABC = arg; };
	inline void setMaxUnitCellEdgeABC(const Double& arg) { MaxUnitCellEdgeABC = arg; };
	inline void setMaxNumPeakInRange(const Int4& arg) { MaxNumPeakInRange = arg; };
	inline void setMinNumPeakInRange(const Int4& arg) { MinNumPeakInRange = arg; };

	inline void setMinFOM(const Double& arg) { MinFOM = arg; };
	inline void setResolution(const Double& arg) { Resol = arg; };
	inline void setMinLatticePointDistance(const Double& arg) { MinLatticePointDistance = arg; };

	// Put functions.
	inline const bool& IsAngleDispersion() const { return IsADorTOF; };
	inline ePeakShiftFunctionType putPeakShiftFunctionType() const;

	inline const bool& putOutputSymmetry(const eBravaisType& i) const { return OutputSymmetry[(Int4)i]; };

	inline const Vec_DP& putConvParam() const{ return ConvParam; };
	inline const Vec_DP& putPeakShiftParamDegree() const{ return PeakShiftParam; };
	inline Vec_DP putPeakShiftParamRadian() const;
	inline const Double& putWaveLength() const{ return WlengthX; };

    // Parameters for search.
	inline const Int4& putMaxPeakNum() const { return MaxPeakNum; };
	inline const Int4& putMaxEdgeNum() const { return MaxEdgeNum; };
	inline const Int4& putMaxNodeNum() const { return MaxNodeNum; };

	inline const Double& putCriticalValueForLinearSum() const { return CVForLinearSum; };
	inline Double putCriticalValueSquareForLinearSum() const { return CVForLinearSum*CVForLinearSum; };
	inline const Double& putMinUnitCellVolume() const { return MinUnitCellVolume; };
	inline const Double& putMaxUnitCellVolume() const { return MaxUnitCellVolume; };

	inline const string& putStrMaxPeakNum() const { return str_MaxPeakNum; };
	inline const string& putStrMaxEdgeNum() const { return str_MaxEdgeNum; };
	inline const string& putStrMinUnitCellVolume() const { return str_MinUnitCellVolume; };
	inline const string& putStrMaxUnitCellVolume() const { return str_MinUnitCellVolume; };
	inline const string& putStrMaxNodeNum() const { return str_MaxNodeNum; };

    // Parameters for output.
	inline const Int4& putNumberOfReflectionsForFigureOfMerit() const { return NumRefFigureOfMerit; };
	inline const Int4& putMaxNumPeakInRange() const { return MaxNumPeakInRange; };
	inline const Int4& putMinNumPeakInRange() const { return MinNumPeakInRange; };

	inline const Double& putMinUnitCellEdgeABC() const { return MinUnitCellEdgeABC; };
	inline const Double& putMaxUnitCellEdgeABC() const { return MaxUnitCellEdgeABC; };
	inline const Double& putMinFOM() const { return MinFOM; };
	inline const Double& putResolution() const { return Resol; };

	inline const string& putStrMaxNumPeakInRange() const { return str_MaxNumPeakInRange; };
	inline const string& putStrMinNumPeakInRange() const { return str_MinNumPeakInRange; };

	inline const eConographAnalysisMode putSearchLevel() const { return (eConographAnalysisMode)SearchLevel; };
	inline const Int4& putNumberOfThreadsToUse() const { return NumCores; };
	inline const Int4& putMaxSizeForEachBRAVAIS() const { return MaxSizeForEachBRAVAIS; };
	inline eRHaxis putRhombohedralAxis() const { return find_key(RHaxisLabel, RhomAxis); };
	inline eABCaxis putBaseCenteredAxis() const { return find_key(ABCaxisLabel, MonoBaseAxis); };
	inline const Double& putThresholdOnNormM() const { return ThresholdNormM; };
	inline const Double& putThresholdOnRevM() const { return ThresholdRevM; };
	inline const Double& putMinLatticePointDistance() const { return MinLatticePointDistance; };
	inline const Double& putResolutionTocheckComputedLines() const { return ResoltoCheckComputedLines; };
	inline const Double& putMaxQTocheckComputedLines() const { return MaxQtoCheckComputedLines; };

	static const string& putMaxPeakNumLabel() { return str_MaxPeakNum_Data.first.putLabel(); };
	static const string& putMinUnitCellVolumeLabel() { return str_MinUnitCellVolume_Data.first.putLabel(); };
	static const string& putMaxUnitCellVolumeLabel() { return str_MaxUnitCellVolume_Data.first.putLabel(); };
	static const string& putNumberOfReflectionsForFigureOfMeritLabel() { return NumRefFigureOfMerit_Data.first.putLabel(); };
	inline const string& putResolutionTocheckComputedLinesLabel() const { return ResoltoCheckComputedLines_Data.first.putLabel(); };
	inline const string& putMaxQTocheckComputedLinesLabel() const { return MaxQtoCheckComputedLines_Data.first.putLabel(); };

	const string& putTagLabel() const;
    void setData(const RWParamProperty& parent_prop,
			vector<RWParam_void>& tray);

	ZErrorMessageReadingFile readFile(const string& filename, const string& file_label);
    ZErrorMessage setAutomaticallyComputedParam(const vector<QData>& qData);
};

inline ePeakShiftFunctionType ControlParam::putPeakShiftFunctionType() const
{
	if( IsADorTOF ) return kPeakShiftFunction_Type1;
	else return kPeakShiftFunction_Type0;
};


inline Vec_DP ControlParam::putPeakShiftParamRadian() const
{
	static const Double RadDeg = PI() / 180.0; // = pi / 180.0.
	Vec_DP ans = PeakShiftParam;
	for(Vec_DP::iterator it=ans.begin(); it<ans.end(); it++) (*it) *= RadDeg;
	return ans;
}

#endif

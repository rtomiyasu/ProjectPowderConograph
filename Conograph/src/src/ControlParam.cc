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
#ifdef _OPENMP
	#include <omp.h>
#endif
#include <string>
#include "ControlParam.hh"
#include "utility_data_structure/qdata.hh"
#include "utility_func/zstring.hh"
#include "utility_func/zmath.hh"
#include "zlog/zlog.hh"

const map<eABCaxis, string> ControlParam::ABCaxisLabel = putABCaxisLabel();
const map<eRHaxis, string> ControlParam::RHaxisLabel = putRHaxisLabel();

static Vec_DP putConvParamRange()
{
	Vec_DP ans(2, 0.0);
	ans[0] = -I_ReadData::MAX_DP();
	return ans;
}

const pair<RWParamProperty, RWParamData<bool> > ControlParam::IsAngleDispersion_Data(
		RWParamProperty(BOOLFLAG, "IsAngleDispersion"),
		RWParamData<bool>(true, REPLACE_NONE<bool>, NULL, false, NULL, false, -1, -1) );
const pair<RWParamProperty, RWParamData<bool> > ControlParam::OutputSymmetry_Data[NUM_LS]
	= {
			pair<RWParamProperty, RWParamData<bool> >( 
				RWParamProperty(BOOLFLAG, "OutputTriclinic"),
				RWParamData<bool>(true, REPLACE_NONE<bool>, NULL, false, NULL, false, -1, -1) ),
			pair<RWParamProperty, RWParamData<bool> >( 
				RWParamProperty(BOOLFLAG, "OutputMonoclinicP"),
				RWParamData<bool>(true, REPLACE_NONE<bool>, NULL, false, NULL, false, -1, -1) ),
			pair<RWParamProperty, RWParamData<bool> >( 
				RWParamProperty(BOOLFLAG, "OutputMonoclinicB"),
				RWParamData<bool>(true, REPLACE_NONE<bool>, NULL, false, NULL, false, -1, -1) ),
			pair<RWParamProperty, RWParamData<bool> >( 
				RWParamProperty(BOOLFLAG, "OutputOrthorhombicP"),
				RWParamData<bool>(true, REPLACE_NONE<bool>, NULL, false, NULL, false, -1, -1) ),
			pair<RWParamProperty, RWParamData<bool> >( 
				RWParamProperty(BOOLFLAG, "OutputOrthorhombicC"),
				RWParamData<bool>(true, REPLACE_NONE<bool>, NULL, false, NULL, false, -1, -1) ),
			pair<RWParamProperty, RWParamData<bool> >( 
				RWParamProperty(BOOLFLAG, "OutputOrthorhombicI"),
				RWParamData<bool>(true, REPLACE_NONE<bool>, NULL, false, NULL, false, -1, -1) ),
			pair<RWParamProperty, RWParamData<bool> >( 
				RWParamProperty(BOOLFLAG, "OutputOrthorhombicF"),
				RWParamData<bool>(true, REPLACE_NONE<bool>, NULL, false, NULL, false, -1, -1) ),
			pair<RWParamProperty, RWParamData<bool> >( 
				RWParamProperty(BOOLFLAG, "OutputTetragonalP"),
				RWParamData<bool>(true, REPLACE_NONE<bool>, NULL, false, NULL, false, -1, -1) ),
			pair<RWParamProperty, RWParamData<bool> >( 
				RWParamProperty(BOOLFLAG, "OutputTetragonalI"),
				RWParamData<bool>(true, REPLACE_NONE<bool>, NULL, false, NULL, false, -1, -1) ),
			pair<RWParamProperty, RWParamData<bool> >( 
				RWParamProperty(BOOLFLAG, "OutputRhombohedral"),
				RWParamData<bool>(true, REPLACE_NONE<bool>, NULL, false, NULL, false, -1, -1) ),
			pair<RWParamProperty, RWParamData<bool> >( 
				RWParamProperty(BOOLFLAG, "OutputHexagonal"),
				RWParamData<bool>(true, REPLACE_NONE<bool>, NULL, false, NULL, false, -1, -1) ),
			pair<RWParamProperty, RWParamData<bool> >( 
				RWParamProperty(BOOLFLAG, "OutputCubicP"),
				RWParamData<bool>(true, REPLACE_NONE<bool>, NULL, false, NULL, false, -1, -1) ),
			pair<RWParamProperty, RWParamData<bool> >( 
				RWParamProperty(BOOLFLAG, "OutputCubicI"),
				RWParamData<bool>(true, REPLACE_NONE<bool>, NULL, false, NULL, false, -1, -1) ),
			pair<RWParamProperty, RWParamData<bool> >( 
				RWParamProperty(BOOLFLAG, "OutputCubicF"),
				RWParamData<bool>(true, REPLACE_NONE<bool>, NULL, false, NULL, false, -1, -1) )
		};

inline Vec_DP put_initial_conv_param()
{
	Vec_DP ans(3,0.0);
	ans[1] = 1.0;
	return ans;
}

const pair<RWParamProperty, RWParamData<Vec_DP> > ControlParam::ConvParam_Data(
		RWParamProperty(DARRAY, "ConversionParameters"), 
		RWParamData<Vec_DP>(put_initial_conv_param(), REPLACE_VECTOR_NONE<Double>, GTVec<Double>, putConvParamRange(), NULL, Vec_DP(), 2, 6) );
const pair<RWParamProperty, RWParamData<Vec_DP> > ControlParam::PeakShiftParam_Data(
		RWParamProperty(DARRAY, "ZeroPointShiftParameter"),
		RWParamData<Vec_DP>(Vec_DP(1,0.0), REPLACE_VECTOR_NONE<Double>, NULL, Vec_DP(), NULL, Vec_DP(), 1, 1) );
const pair<RWParamProperty, RWParamData<Double> > ControlParam::WlengthX_Data(
		RWParamProperty(DVALUE, "WaveLength"), 
		RWParamData<Double>(1.54056, REPLACE_NONE<Double>, GT<Double>, 0.0, NULL, MAX_DP(), -1, -1) );


const pair<RWParamProperty, RWParamData<string> > ControlParam::str_MaxPeakNum_Data(
		RWParamProperty(STRVALUE, "MaxNumberOfPeaks"),
		RWParamData<string>("AUTO", REPLACE_NONE<string>, NULL, "", NULL, "", -1, -1) );
const pair<RWParamProperty, RWParamData<Int4> > ControlParam::NumRefFigureOfMerit_Data(
		RWParamProperty(INT4VALUE, "MaxNumberOfPeaksForFOM"), 
		RWParamData<Int4>(20, REPLACE_NONE<Int4>, GE<Int4>, 1, NULL, MAX_INT(), -1, -1) );	// 20 <= param < INF.
const pair<RWParamProperty, RWParamData<string> > ControlParam::str_MaxEdgeNum_Data(
		RWParamProperty(STRVALUE, "MaxNumberOfTwoDimTopographs"),
		RWParamData<string>("AUTO", REPLACE_NONE<string>, NULL, "", NULL, "", -1, -1) );	// 0 < param < INF.
const pair<RWParamProperty, RWParamData<string> > ControlParam::str_MaxNodeNum_Data(
		RWParamProperty(STRVALUE, "MaxNumberOfLatticeCandidates"),
		RWParamData<string>("AUTO", REPLACE_NONE<string>, NULL, "", NULL, "", -1, -1) );	// 0 < param < INF.
const pair<RWParamProperty, RWParamData<string> > ControlParam::str_MaxNumPeakInRange_Data(
		RWParamProperty(STRVALUE, "MaxNumberOfMillerIndicesInRange"),
		RWParamData<string>("AUTO", REPLACE_NONE<string>, NULL, "", NULL, "", -1, -1) );
const pair<RWParamProperty, RWParamData<string> > ControlParam::str_MinNumPeakInRange_Data(
		RWParamProperty(STRVALUE, "MinNumberOfMillerIndicesInRange"),
		RWParamData<string>("AUTO", REPLACE_NONE<string>, NULL, "", NULL, "", -1, -1) );	// 0.0 <= param <= 1.0.
const pair<RWParamProperty, RWParamData<Double> > ControlParam::MinFOM_Data(
		RWParamProperty(DVALUE, "MinFOM"), 
		RWParamData<Double>(3.0, REPLACE_NONE<Double>, GE<Double>, 0.0, NULL, MAX_DP(), -1, -1) );	// 0 <= param < INF.
const pair<RWParamProperty, RWParamData<string> > ControlParam::str_MaxUnitCellVolume_Data(
		RWParamProperty(STRVALUE, "MaxPrimitiveUnitCellVolume"),
		RWParamData<string>("AUTO", REPLACE_NONE<string>, NULL, "", NULL, "", -1, -1) );	// 0 <= param < INF.
const pair<RWParamProperty, RWParamData<string> > ControlParam::str_MinUnitCellVolume_Data(
		RWParamProperty(STRVALUE, "MinPrimitiveUnitCellVolume"),
		RWParamData<string>("AUTO", REPLACE_NONE<string>, NULL, "", NULL, "", -1, -1) );	// 0 <= param < INF.
const pair<RWParamProperty, RWParamData<Double> > ControlParam::MinUnitCellEdgeABC_Data(
		RWParamProperty(DVALUE, "MinUnitCellEdgeABC"),
		RWParamData<Double>(0.0, REPLACE_NONE<Double>, GE<Double>, 0.0, NULL, MAX_DP(), -1, -1) );	// 0 <= param < INF.
const pair<RWParamProperty, RWParamData<Double> > ControlParam::MaxUnitCellEdgeABC_Data(
		RWParamProperty(DVALUE, "MaxUnitCellEdgeABC"),
		RWParamData<Double>(1000.0, REPLACE_MAX<Double>, GT<Double>, 0.0, NULL, MAX_DP(), -1, -1) );	// 0 <= param < INF.
const pair<RWParamProperty, RWParamData<Double> > ControlParam::CVForLinearSum_Data(
		RWParamProperty(DVALUE, "CriticalValueForLinearSum"), 
		RWParamData<Double>(1.0, REPLACE_NONE<Double>, GT<Double>, 0.0, NULL, MAX_DP(), -1, -1) );	// 0 < param < INF.
const pair<RWParamProperty, RWParamData<Double> > ControlParam::Resol_Data(
		RWParamProperty(DVALUE, "Resolution"), 
		RWParamData<Double>(0.03, REPLACE_NONE<Double>, GE<Double>, 0.0, LE<Double>, 0.25, -1, -1) );	// 0 <= param < INF.

const pair<RWParamProperty, RWParamData<Int4> > ControlParam::SearchLevel_Data(
		RWParamProperty(INT4VALUE, "SearchLevel"),
		RWParamData<Int4>((Int4)ConographQuickSearch, REPLACE_NONE<Int4>, GE<Int4>, (Int4)ConographQuickSearch, LE<Int4>, (Int4)ConographRegularSearch, -1, -1) );

const pair<RWParamProperty, RWParamData<string> > ControlParam::MonoBaseAxis_Data(
		RWParamProperty(STRVALUE, "AxisForMonoclinicSymmetry"),
		RWParamData<string>("B", REPLACE_NONE<string>, NULL, "", NULL, "", -1, -1) );
const pair<RWParamProperty, RWParamData<string> > ControlParam::RhomAxis_Data(
		RWParamProperty(STRVALUE, "AxisForRhombohedralSymmetry"),
		RWParamData<string>("Hexagonal", REPLACE_NONE<string>, NULL, "", NULL, "", -1, -1) );
const pair<RWParamProperty, RWParamData<Double> > ControlParam::ThresholdNormM_Data(
		RWParamProperty(DVALUE, "ThresholdOnNormM"),
		RWParamData<Double>(1.9, REPLACE_NONE<Double>, GE<Double>, 0.0, NULL, MAX_DP(), -1, -1) );	// 0 <= param < INF.
const pair<RWParamProperty, RWParamData<Double> > ControlParam::ThresholdRevM_Data(
		RWParamProperty(DVALUE, "ThresholdOnRevM"),
		RWParamData<Double>(1.0, REPLACE_NONE<Double>, GE<Double>, 0.0, NULL, MAX_DP(), -1, -1) );	// 0 <= param < INF.
const pair<RWParamProperty, RWParamData<Double> > ControlParam::MinLatticePointDistance_Data(
		RWParamProperty(DVALUE, "MinDistanceBetweenLatticePoints"),
		RWParamData<Double>(2.0, REPLACE_NONE<Double>, GE<Double>, 0.0, NULL, MAX_DP(), -1, -1) );	// 0 <= param < INF.
const pair<RWParamProperty, RWParamData<Double> > ControlParam::ResoltoCheckComputedLines_Data(
		RWParamProperty(DVALUE, "ResolutionToCheckCoincidenceofComputedLines", NULL, 1, 0),
		RWParamData<Double>(0.03, REPLACE_NONE<Double>, GE<Double>, 0.0, LE<Double>, 0.05, -1, -1) );	// 0 <= param < INF.
const pair<RWParamProperty, RWParamData<Double> > ControlParam::MaxQtoCheckComputedLines_Data(
		RWParamProperty(DVALUE, "MaxQToCheckCoincidenceofComputedLines", NULL, 1, 0),
		RWParamData<Double>(3.0, REPLACE_NONE<Double>, GE<Double>, 0.0, NULL, MAX_DP(), -1, -1) );	// 0 <= param < INF.
const pair<RWParamProperty, RWParamData<Int4> > ControlParam::MaxSizeForEachBRAVAIS_Data(
		RWParamProperty(INT4VALUE, "MaxNumberOfSolutionsForEachBravaisLattice", NULL, 1, 0),
		RWParamData<Int4>(1000, REPLACE_NONE<Int4>, GE<Int4>, 1, NULL, MAX_INT(), -1, -1) );

ControlParam::ControlParam()
	:
#ifdef _OPENMP
		NumCores_Data(RWParamProperty(INT4VALUE, "NumberOfThreadsToUse"),
		RWParamData<Int4>(max(1, omp_get_num_procs()-1), REPLACE_NUM_THREAD, GE<Int4>, 1, NULL, MAX_INT(), -1, -1) ),
#else
		NumCores_Data(RWParamProperty(INT4VALUE, "NumberOfThreadsToUse"),
		RWParamData<Int4>(1, REPLACE_NUM_THREAD, GE<Int4>, 1, NULL, MAX_INT(), -1, -1) ),
#endif
		ReadConfigurationParameters(true),
		IsADorTOF(IsAngleDispersion_Data.second.initial_value),
		ConvParam(ConvParam_Data.second.initial_value),
		PeakShiftParam(PeakShiftParam_Data.second.initial_value),
		WlengthX(WlengthX_Data.second.initial_value),
		SearchLevel(SearchLevel_Data.second.initial_value),
		MaxPeakNum(0),
		MaxEdgeNum(0),
		MaxNodeNum(0),
		MinUnitCellVolume(0.0),
		MaxUnitCellVolume(0.0),
		CVForLinearSum(CVForLinearSum_Data.second.initial_value),
		NumRefFigureOfMerit(NumRefFigureOfMerit_Data.second.initial_value),
		MaxNumPeakInRange(0),
		MinNumPeakInRange(0),
		MinUnitCellEdgeABC(MinUnitCellEdgeABC_Data.second.initial_value),
		MaxUnitCellEdgeABC(MaxUnitCellEdgeABC_Data.second.initial_value),
		MinFOM(MinFOM_Data.second.initial_value),
		Resol(Resol_Data.second.initial_value),
		NumCores(NumCores_Data.second.initial_value),
		MaxSizeForEachBRAVAIS(MaxSizeForEachBRAVAIS_Data.second.initial_value),
		MonoBaseAxis(MonoBaseAxis_Data.second.initial_value),
		RhomAxis(RhomAxis_Data.second.initial_value),
		ThresholdNormM(ThresholdNormM_Data.second.initial_value),
		ThresholdRevM(ThresholdRevM_Data.second.initial_value),
		MinLatticePointDistance(MinLatticePointDistance_Data.second.initial_value),
		ResoltoCheckComputedLines(ResoltoCheckComputedLines_Data.second.initial_value),
		MaxQtoCheckComputedLines(MaxQtoCheckComputedLines_Data.second.initial_value),
		str_MaxPeakNum(str_MaxPeakNum_Data.second.initial_value),
		str_MaxEdgeNum(str_MaxEdgeNum_Data.second.initial_value),
		str_MaxNodeNum(str_MaxNodeNum_Data.second.initial_value),
		str_MinUnitCellVolume(str_MinUnitCellVolume_Data.second.initial_value),
		str_MaxUnitCellVolume(str_MaxUnitCellVolume_Data.second.initial_value),
		str_MaxNumPeakInRange(str_MaxNumPeakInRange_Data.second.initial_value),
		str_MinNumPeakInRange(str_MinNumPeakInRange_Data.second.initial_value)
{
	for(Int4 i=0; i<NUM_LS; i++)
	{
		OutputSymmetry[i] = OutputSymmetry_Data[i].second.initial_value;
	}
}


ControlParam::~ControlParam()
{
}


const string& ControlParam::putTagLabel() const
{
	static const string label = "ConographParameters";
	return label;
}


ZErrorMessageReadingFile ControlParam::readFile(const string& filename, const string& file_label)
{
	ZErrorMessageReadingFile zerr = I_ReadData::readFile(filename, file_label);
	if( IsAngleDispersion() ) ConvParam.clear();
	else PeakShiftParam.clear();

	if( zerr.putErrorType() != ZErrorNoError ) return zerr;
	if( str_MaxPeakNum!="AUTO" )
	{
		if( str_MaxPeakNum=="MAX" )
		{
			MaxPeakNum = numeric_limits<Int4>::max();
		}
		else
		{
			istringstream iss( str_MaxPeakNum );
			iss >> MaxPeakNum;
			if( iss.fail() ) return ZErrorMessageReadingFile( filename, ZErrorMessage(ZErrorFileFormatBroken, "The string cannot be replaced for <" + str_MaxPeakNum_Data.first.putLabel() + ">:"+str_MaxPeakNum, __FILE__, __LINE__, __FUNCTION__) );
			if( MaxPeakNum < 3 )
			{
				return ZErrorMessageReadingFile( filename, ZErrorMessage("<"+ str_MaxPeakNum_Data.first.putLabel() + "> is too small: "+num2str<Int4>(MaxPeakNum), __FILE__, __LINE__, __FUNCTION__) );
			}
			else if( MaxPeakNum > 2000 )
			{
				return ZErrorMessageReadingFile( filename, ZErrorMessage("<"+ str_MaxPeakNum_Data.first.putLabel() + "> is too large: "+num2str<Int4>(MaxPeakNum), __FILE__, __LINE__, __FUNCTION__) );
			}
		}
	}

	if( str_MaxEdgeNum != "AUTO" )
	{
		if( str_MaxEdgeNum=="MAX" )
		{
			MaxEdgeNum = numeric_limits<Int4>::max();
		}
		else
		{
			istringstream iss( str_MaxEdgeNum );
			iss >> MaxEdgeNum;
			if( iss.fail() ) return ZErrorMessageReadingFile(filename, ZErrorMessage(ZErrorFileFormatBroken, "The string cannot be replaced for <" + str_MaxEdgeNum_Data.first.putLabel() + ">: "+str_MaxEdgeNum, __FILE__, __LINE__, __FUNCTION__) );
			if( MaxEdgeNum < 1 )
			{
				return ZErrorMessageReadingFile(filename, nerror_arg("<"+ str_MaxEdgeNum_Data.first.putLabel() + "> is too small: "+num2str<Int4>(MaxEdgeNum), __FILE__, __LINE__, __FUNCTION__) );
			}
		}
	}

	if( str_MaxNodeNum != "AUTO" )
	{
		if( str_MaxNodeNum=="MAX" )
		{
			MaxNodeNum = numeric_limits<Int4>::max();
		}
		else
		{
			istringstream iss( str_MaxNodeNum );
			iss >> MaxNodeNum;
			if( iss.fail() ) return ZErrorMessageReadingFile(filename, ZErrorMessage(ZErrorFileFormatBroken, "The string in <" + str_MaxNodeNum_Data.first.putLabel() + "> is not replaced: "+str_MaxNodeNum, __FILE__, __LINE__, __FUNCTION__) );
			if( MaxNodeNum < 1 )
			{
				return ZErrorMessageReadingFile(filename, nerror_arg("<"+ str_MaxNodeNum_Data.first.putLabel() + "> is too small: "+num2str<Int4>(MaxNodeNum), __FILE__, __LINE__, __FUNCTION__) );
			}
		}
	}

	if( str_MinUnitCellVolume != "AUTO" )
	{
		istringstream iss( str_MinUnitCellVolume );
		iss >> MinUnitCellVolume;
		if( iss.fail() ) return ZErrorMessageReadingFile(filename, ZErrorMessage(ZErrorFileFormatBroken, "The string cannot be replaced for <" + putMinUnitCellVolumeLabel() + ">: "+str_MinUnitCellVolume, __FILE__, __LINE__, __FUNCTION__) );
		if( MinUnitCellVolume < 1.0 )
		{
			return ZErrorMessageReadingFile(filename, nerror_arg("<"+ putMinUnitCellVolumeLabel() + "> is less than 1: "+num2str<Double>(MinUnitCellVolume), __FILE__, __LINE__, __FUNCTION__) );
		}
	}

	if( str_MaxUnitCellVolume != "AUTO" )
	{
		istringstream iss( str_MaxUnitCellVolume );
		iss >> MaxUnitCellVolume;
		if( iss.fail() ) return ZErrorMessageReadingFile(filename, ZErrorMessage(ZErrorFileFormatBroken, "The string cannot be replaced for <" + putMaxUnitCellVolumeLabel() + ">: "+str_MaxUnitCellVolume, __FILE__, __LINE__, __FUNCTION__) );
		if( MaxUnitCellVolume <= MinUnitCellVolume )
		{
			return ZErrorMessageReadingFile(filename, nerror_arg("<"+ putMaxUnitCellVolumeLabel() + "> is not larger than <"+ putMinUnitCellVolumeLabel() + "> : "+num2str<Double>(MaxUnitCellVolume), __FILE__, __LINE__, __FUNCTION__) );
		}
	}

	if( str_MaxNumPeakInRange != "AUTO" )
	{
		if( str_MaxNumPeakInRange =="MAX" )
		{
			MaxNumPeakInRange = numeric_limits<Int4>::max();
		}
		else
		{
			istringstream iss( str_MaxNumPeakInRange );
			iss >> MaxNumPeakInRange;
			if( iss.fail() ) return ZErrorMessageReadingFile(filename, ZErrorMessage(ZErrorFileFormatBroken, "The string cannot be replaced for <" + str_MaxNumPeakInRange_Data.first.putLabel() + ">: "+str_MaxNumPeakInRange, __FILE__, __LINE__, __FUNCTION__) );
			if( MaxNumPeakInRange < 1 )
			{
				return ZErrorMessageReadingFile(filename, ZErrorMessage("<"+ str_MaxNumPeakInRange_Data.first.putLabel() + "> is less than 1: "+num2str<Double>(MaxNumPeakInRange), __FILE__, __LINE__, __FUNCTION__) );
			}
		}
	}

	if( str_MinNumPeakInRange != "AUTO" )
	{
		istringstream iss( str_MinNumPeakInRange );
		iss >> MinNumPeakInRange;
		if( iss.fail() ) return ZErrorMessageReadingFile(filename, ZErrorMessage(ZErrorFileFormatBroken, "The string cannot be replaced for <" + str_MinNumPeakInRange_Data.first.putLabel() + ">: "+str_MinNumPeakInRange, __FILE__, __LINE__, __FUNCTION__) );
		if( MinNumPeakInRange < 1 )
		{
			return ZErrorMessageReadingFile(filename, ZErrorMessage("<"+ str_MinNumPeakInRange_Data.first.putLabel() + "> is less than 1: "+num2str<Double>(MinNumPeakInRange), __FILE__, __LINE__, __FUNCTION__) );
		}

		if( str_MaxNumPeakInRange != "AUTO" && MaxNumPeakInRange < MinNumPeakInRange )
		{
			return ZErrorMessageReadingFile(filename, ZErrorMessage( "<" + str_MaxNumPeakInRange_Data.first.putLabel() + "> is less than <" + str_MinNumPeakInRange_Data.first.putLabel() + ">",
											__FILE__, __LINE__, __FUNCTION__) );
		}
	}

	return ZErrorMessageReadingFile();
}



void ControlParam::setData(const RWParamProperty& parent_prop,
		vector<RWParam_void>& tray)
{
	if( parent_prop.putLabel() == this->putTagLabel() )
	{
		tray.push_back( RWParam_void(IsAngleDispersion_Data, &IsADorTOF) );
		for(Int4 i=0; i<NUM_LS; i++)
		{
			tray.push_back( RWParam_void(OutputSymmetry_Data[i], &OutputSymmetry[i]) );
		}
		tray.push_back( RWParam_void(ConvParam_Data, &ConvParam) );
		tray.push_back( RWParam_void(PeakShiftParam_Data, &PeakShiftParam) );
		tray.push_back( RWParam_void(WlengthX_Data, &WlengthX) );
		tray.push_back( RWParam_void(Resol_Data, &Resol) );
		tray.push_back( RWParam_void(CVForLinearSum_Data, &CVForLinearSum) );
		tray.push_back( RWParam_void(str_MaxPeakNum_Data, &str_MaxPeakNum) );
		tray.push_back( RWParam_void(str_MaxEdgeNum_Data, &str_MaxEdgeNum) );
		tray.push_back( RWParam_void(str_MaxNodeNum_Data, &str_MaxNodeNum) );
		tray.push_back( RWParam_void(str_MinUnitCellVolume_Data, &str_MinUnitCellVolume) );
		tray.push_back( RWParam_void(str_MaxUnitCellVolume_Data, &str_MaxUnitCellVolume) );

		tray.push_back( RWParam_void(MinUnitCellEdgeABC_Data, &MinUnitCellEdgeABC) );
		tray.push_back( RWParam_void(MaxUnitCellEdgeABC_Data, &MaxUnitCellEdgeABC) );
		tray.push_back( RWParam_void(MinFOM_Data, &MinFOM) );
		tray.push_back( RWParam_void(str_MaxNumPeakInRange_Data, &str_MaxNumPeakInRange) );
		tray.push_back( RWParam_void(str_MinNumPeakInRange_Data, &str_MinNumPeakInRange) );
		tray.push_back( RWParam_void(NumRefFigureOfMerit_Data, &NumRefFigureOfMerit) );

		tray.push_back( RWParam_void(SearchLevel_Data, &SearchLevel) );

		if( ReadConfigurationParameters )
		{
			tray.push_back( RWParam_void(NumCores_Data, &NumCores) );
			MaxSizeForEachBRAVAIS = MaxSizeForEachBRAVAIS_Data.second.initial_value;
			tray.push_back( RWParam_void(MaxSizeForEachBRAVAIS_Data, &MaxSizeForEachBRAVAIS) );

			tray.push_back( RWParam_void(MonoBaseAxis_Data, &MonoBaseAxis) );
			tray.push_back( RWParam_void(RhomAxis_Data, &RhomAxis) );
			tray.push_back( RWParam_void(ThresholdNormM_Data, &ThresholdNormM) );
			tray.push_back( RWParam_void(ThresholdRevM_Data, &ThresholdRevM) );
			tray.push_back( RWParam_void(MinLatticePointDistance_Data, &MinLatticePointDistance) );
			tray.push_back( RWParam_void(ResoltoCheckComputedLines_Data, &ResoltoCheckComputedLines) );
			tray.push_back( RWParam_void(MaxQtoCheckComputedLines_Data, &MaxQtoCheckComputedLines) );
		}
	}
}


ZErrorMessage ControlParam::checkData(const RWParam_void& param) const
{
	const string Label = param.putLabel(); 
	if( IsEqualTag(param.putProperty(), ConvParam_Data.first ) )
	{
		if( IsADorTOF ) 
		{
			return ZErrorMessage();
		}
	}
	else if( IsEqualTag(param.putProperty(), PeakShiftParam_Data.first)
			|| IsEqualTag(param.putProperty(), WlengthX_Data.first) )
	{
		if( !IsADorTOF ) 
		{
			return ZErrorMessage();
		}
	}
	else if( IsEqualTag(param.putProperty(), MonoBaseAxis_Data.first) )
	{
		if( find_key(ABCaxisLabel, MonoBaseAxis) == eABCaxis(-1) )
		{
			return nerror_out_range(param.putLabel(), MonoBaseAxis, __FILE__, __LINE__, __FUNCTION__ );
		}
		else return ZErrorMessage();
	}
	else if( IsEqualTag(param.putProperty(), RhomAxis_Data.first) )
	{
		if( find_key(RHaxisLabel, RhomAxis) == eRHaxis(-1) )
		{
			return nerror_out_range(param.putLabel(), RhomAxis, __FILE__, __LINE__, __FUNCTION__ );
		}
		else return ZErrorMessage();
	}
	return I_ReadData::checkData(param);
}



ZErrorMessage ControlParam::checkIfDataAreSet(const RWParam_void& param, const Int4& num) const
{
	const string Label = param.putLabel();
	if( IsEqualTag(param.putProperty(), ConvParam_Data.first) )
	{
		if( IsADorTOF ) 
		{
			return ZErrorMessage();
		}
	}
	else if( IsEqualTag(param.putProperty(), PeakShiftParam_Data.first)
			|| IsEqualTag(param.putProperty(), WlengthX_Data.first) )
	{
		if( !IsADorTOF ) 
		{
			return ZErrorMessage();
		}
	}
	else if( IsEqualTag(param.putProperty(), MaxUnitCellEdgeABC_Data.first) )
	{
		if( MaxUnitCellEdgeABC < MinUnitCellEdgeABC )
		{
			return ZErrorMessage( "<" + MaxUnitCellEdgeABC_Data.first.putLabel() + "> is less than <" + MinUnitCellEdgeABC_Data.first.putLabel() + ">",
										__FILE__, __LINE__, __FUNCTION__);
		}
	}

	return I_ReadData::checkIfDataAreSet(param, num);
}


inline Int4 put_automatic_peak_num(const eConographAnalysisMode& search_level, const vector<QData>& qdata,
		const Int4& num_ref_figure_of_merit, const Double& min_length_edge)
{
	assert( num_ref_figure_of_merit <= (Int4)qdata.size() );
	if( qdata.empty() ) return 0;
	return max(num_ref_figure_of_merit, min(48, (Int4)distance( qdata.begin(),
			upper_bound( qdata.begin(), qdata.end(), QData( 10.0 /(min_length_edge*min_length_edge), 0.0 ) ) ) ) );
}


inline Int4 put_automatic_edge_num(const eConographAnalysisMode& search_level, const Int4& MaxPeakNum)
{
	if( search_level == ConographQuickSearch ) return max(100, MaxPeakNum*2);
	return MaxPeakNum * (MaxPeakNum+1) / 2;
}


inline Int4 put_automatic_node_num(const eConographAnalysisMode& search_level, const Int4& MaxPeakNum)
{
	if( search_level == ConographQuickSearch ) return 64000;
	return min(MaxPeakNum * (MaxPeakNum+1) * (MaxPeakNum+2) * 2 / 3, 32000);
}


inline pair<Double, Double> put_automatic_unit_cell_volume(const vector<QData>& qData, const Int4& MaxPeakNum)
{
	Double ans = 1.0/put_maximum_counitcell_volume(qData, 20);
	return pair<Double, Double>(max(5.0, ans), ans*30.0);
}


ZErrorMessage ControlParam::setAutomaticallyComputedParam(const vector<QData>& qData)
{
	// Set the automatically calculated number of peaks to utilize for indexing.
	if( this->putNumberOfReflectionsForFigureOfMerit() > (Int4)qData.size() )
	{
		this->setNumberOfReflectionsForFigureOfMerit(qData.size());
ZLOG_INFO( "\"" +  ControlParam::putNumberOfReflectionsForFigureOfMeritLabel() + "\" is more than the number of peaks.  It is replaced by the number of peaks.\n" );
	}

	if( this->putStrMaxPeakNum()=="AUTO" )
	{
		MaxPeakNum = put_automatic_peak_num( this->putSearchLevel(), qData, NumRefFigureOfMerit, MinLatticePointDistance);
ZLOG_INFO( "<" + str_MaxPeakNum_Data.first.putLabel() + "> is set to " + num2str<Int4>(MaxPeakNum) + "\n" );
	}
	else
	{
		if( this->putMaxPeakNum() > (Int4)qData.size() )
		{
			this->setMaxPeakNum(qData.size());
ZLOG_INFO( "\"" +  str_MaxPeakNum_Data.first.putLabel() + "\" is more than the number of peaks.  It is replaced by the number of peaks.\n" );
		}
	}

	if( str_MaxEdgeNum == "AUTO" )
	{
		MaxEdgeNum = put_automatic_edge_num( this->putSearchLevel(), this->putMaxPeakNum() );
ZLOG_INFO( "<" + str_MaxEdgeNum_Data.first.putLabel() + "> is set to " + num2str(MaxEdgeNum) + ".\n" );
	}

	if( str_MaxNodeNum == "AUTO" )
	{
		MaxNodeNum = put_automatic_node_num( this->putSearchLevel(), this->putMaxPeakNum() );
ZLOG_INFO( "<" + str_MaxNodeNum_Data.first.putLabel() + "> is set to " + num2str(MaxNodeNum) + ".\n" );
	}

	if( str_MinUnitCellVolume == "AUTO" || str_MaxUnitCellVolume == "AUTO" )
	{
		pair<Double, Double> unit_cell_vol_range = put_automatic_unit_cell_volume(qData, this->putMaxPeakNum());
		if( str_MinUnitCellVolume == "AUTO" )
		{
			MinUnitCellVolume = unit_cell_vol_range.first;
ZLOG_INFO( "<" + putMinUnitCellVolumeLabel() + "> is set to " + num2str(MinUnitCellVolume) + ".\n" );
		}
		if( str_MaxUnitCellVolume == "AUTO" )
		{
			MaxUnitCellVolume = unit_cell_vol_range.second;
ZLOG_INFO( "<" + putMaxUnitCellVolumeLabel() + "> is set to " + num2str(MaxUnitCellVolume) + ".\n" );
		}
	}

	if( str_MaxNumPeakInRange == "AUTO" )
	{
		MaxNumPeakInRange = ifloor( pow(this->putNumberOfReflectionsForFigureOfMerit(), 1.62) );
ZLOG_INFO( "<" + str_MaxNumPeakInRange_Data.first.putLabel() + "> is set to " + num2str(MaxNumPeakInRange) + ".\n" );
	}

	if( str_MinNumPeakInRange == "AUTO" )
	{
		MinNumPeakInRange = this->putNumberOfReflectionsForFigureOfMerit() * 3 / 5;
ZLOG_INFO( "<" + str_MinNumPeakInRange_Data.first.putLabel() + "> is set to " + num2str(MinNumPeakInRange) + ".\n" );
	}

	if( MaxNumPeakInRange < MinNumPeakInRange )
	{
		return ZErrorMessage("<" + str_MinNumPeakInRange_Data.first.putLabel() + "> : " + num2str<Int4>(MinNumPeakInRange) +'\n'
								 + str_MaxNumPeakInRange_Data.first.putLabel() + "> : " + num2str<Int4>(MaxNumPeakInRange) +'\n'
								 + "<" + str_MaxNumPeakInRange_Data.first.putLabel() + "> is less than <" + str_MinNumPeakInRange_Data.first.putLabel() + ">", __FILE__, __LINE__, __FUNCTION__);
	}

	return ZErrorMessage();
}

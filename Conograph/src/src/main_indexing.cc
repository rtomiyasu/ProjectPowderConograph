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
#include <iostream>
#include <string>
#include <stdexcept>
#include <cstdlib>
#include <time.h>
#include "ControlFile.hh"
#include "ControlParam.hh"
#include "IndexingLattice.hh"
#include "SortingLattice.hh"
#include "lattice_symmetry/OutputInfo.hh"
#include "lattice_symmetry/VCLatticeFigureOfMeritToCheckSymmetry.hh"
#include "lattice_symmetry/LatticeFigureOfMeritToDisplay.hh"
#include "LatticeWithSameQ/LatticeMetricTensor.hh"
#include "LatticeWithSameQ/LatticeWithSameQ.hh"
#include "LatticeWithSameQ/p_out_same_q.hh"
#include "zerror_type/error_out.hh"
#include "zerror_type/error_mes.hh"
#include "zlog/zlog.hh"
#include "p_out_indexing.hh"
#include "qc/p_out_space_group_dtm.hh"
#include "chToqValue.hh"
#include "utility_func/stopx.hh"
#include "utility_func/zstring.hh"


using namespace std;


int main(int argc, char* argv[])
{
	clock_t start = clock();    /* Record the starting time. */

	// For interruption signal.
	SetSignal(SIGINT);

	static const string controlFile = "cntl.inp.xml";
	static const string InputFileLabel = "ZCodeParameters";
    
	try{
	    CRLog::append(new CCoutListner());
	    CRLog::append(new FileoutListner("LOG_CONOGRAPH.txt", zListnerID(1)));

ZLOG_INFO( "Reading " + controlFile + "...\n\n" );

		ControlFile cf;
    	ZErrorMessage zerr;
    	zerr = cf.readFile(controlFile, InputFileLabel);
    	if( zerr.putErrorType() != ZErrorNoError) throw zerr;

    	PeakPosData pData;
    	zerr = pData.readFile(cf.putPeakDataFileName());
    	if( zerr.putErrorType() != ZErrorNoError) throw zerr;

    	if( pData.putPeakPosXData().size() < 3 )
    	{
    		throw nerror_arg("NUMBER OF INPUT REFLECTIONS IS TOO SMALL.", __FILE__, __LINE__, __FUNCTION__);
    	}

    	ControlParam cData;
    	zerr = cData.readFile(cf.putControlParamFileName(), InputFileLabel);
    	if( zerr.putErrorType() != ZErrorNoError ) throw zerr;

#ifdef _OPENMP
    	omp_set_num_threads(min(omp_get_max_threads(), cData.putNumberOfThreadsToUse()));
ZLOG_INFO( "The number of threads is set to " + num2str( min(omp_get_max_threads(), cData.putNumberOfThreadsToUse()) ) + "\n" );
#endif

		// change the peak-positions into q-values.
    	vector<QData> qData;
    	vector<Int4> qIndexData;
    	zerr = chToqValue(cData, pData, qData, qIndexData);
    	if( zerr.putErrorType() != ZErrorNoError) throw zerr;

    	// Set Qdata as a static member of VCData.
    	VCData::setQData(qData, qIndexData);

#ifdef DEBUG
ZLOG_INFO( "Outputting q-values...\nNo., qvalue, error_of_qvalue, peakpos, peak-width, flag\n" );
	stringstream strstream;
	for(UInt4 k=0; k<qData.size(); k++)
	{
		strstream << qIndexData[k] + 1 << "  "
			 << qData[k].q << "  "
			 << sqrt( qData[k].q_var ) << "  "
			 << (pData.putPeakPosXData())[ qIndexData[k] ] << "  "
			 << (pData.putPeakWidthData())[ qIndexData[k] ] << "  "
			 << (pData.putToUseFlag())[ qIndexData[k] ] << endl;
	}
ZLOG_INFO( strstream.str() + "\n" );
#endif

		if( (Int4)qData.size() < 3 )
		{
			throw ZErrorMessage("The number of q-values is too small : " + num2str<Int4>(qData.size()), __FILE__, __LINE__, __FUNCTION__);
		}

		// Check if the range of q-values is large enough.
		zerr = cData.setAutomaticallyComputedParam(qData);
	   	if( zerr.putErrorType() != ZErrorNoError ) throw zerr;
	   	const Double inv_d_max = sqrt( min( qData.begin() + cData.putMaxPeakNum(), qData.end()-1 )->q );
	   	const Double inv_d_min = sqrt( qData.begin()->q );
ZLOG_INFO( "\nThe range of d^*-value : sqrt(the minimum of q-values)--sqrt(the maximum of q-values) = "
	 + num2str( inv_d_min ) + "--" + num2str( inv_d_max ) + "\n" );

		// Estimation of zero point shift.
		if( cData.IsAngleDispersion()
    		&& cData.putPeakShiftFunctionType() !=  kPeakShiftFunction_Type0 )
		{
    		vector< pair<Int4, Int4> > pair_of_q;
			vector< Double > zero_point_shift_deg;
			fitZeroPointShift(pData, 10, pair_of_q, zero_point_shift_deg);

ZLOG_INFO( "-------- Estimated zero point shift from pairs of peaks --------\nNo. of peaks (ratio of sin(theta)) : zero point shift\n" );
			const Int4 ISIZE = zero_point_shift_deg.size();
			if( ISIZE > 0 )
			{
				const Vec_DP& posdata = pData.putPeakPosXData();
				static const Double RadDeg = PI() / 180.0; // = pi / 180.0.
				stringstream strstream;
				strstream.precision(4);
				strstream.setf(ios::fixed, ios::floatfield);
				for(Int4 i=0; i<ISIZE; i++)
				{
					strstream.width(4);
					strstream << num2str( pair_of_q[i].first + 1 ) + "&";
					strstream.width(4);
					strstream << num2str( pair_of_q[i].second + 1 ) + " ("
						 << sin(0.5*posdata[ pair_of_q[i].second ]*RadDeg) / sin(0.5*posdata[ pair_of_q[i].first ]*RadDeg) << "): ";
					strstream.width(7);
					strstream << zero_point_shift_deg[i] << "\n";
				}
				ZLOG_INFO( strstream.str() );
			}
			else
			{
				ZLOG_INFO( "Could not estimate because the number of peaks is small\n" );
			}
ZLOG_INFO( "-------- Estimated zero point shift from pairs of peaks --------\n" );
		}

    	// Indexing.
    	IndexingLattice clc;
    	clc.setParam(cData);

    	clc.determineTwoDimLattices(pData, "DEBUG_CONOGRAPH");
    	clc.determineThreeDimLattices("DEBUG_CONOGRAPH");
    	const vector<SymMat43_VCData>& S_super = clc.putThreeDimTopographNodes();

ZLOG_INFO( "The program has obtained " + num2str( S_super.size() )
		+ " unit-cell parameters in CPU time : " + num2str( (clock() - start) / CLOCKS_PER_SEC )
		+ " [sec.]\n\n" );

		// Indexing.
		static const Int4 NUM_LS = put_number_of_bravais_types();
		vector<VCLatticeFigureOfMeritToCheckSymmetry> lattice_result[NUM_LS];
		
		SortingLattice srl;
		srl.setParam(cData);

		start = clock();    /* Record the starting time. */
		srl.putLatticeCandidatesForEachBravaisTypes(S_super, cData.putThresholdOnNormM(), cData.putThresholdOnRevM(),
										cData.putMaxSizeForEachBRAVAIS(),
										cData.putBaseCenteredAxis(), cData.putRhombohedralAxis(), lattice_result);

ZLOG_INFO( "Selecting lattices with the best figure of merit among almost equivalent solutions...\n" );
		// After this method, each lattice_result[i] are sorted by unit-cell volumes.
		srl.setNumberOfNeighbors(cData.putBaseCenteredAxis(), OutputInfo::CmpFunc[(Int4)SCM], lattice_result);

ZLOG_INFO( "The Bravais lattice determination and refinement of unit-cell parameters finish in CPU time : " + num2str( (clock() - start) / CLOCKS_PER_SEC ) + " [sec.]\n\n" );

		// Sort by sort_criterion.
		static const eSortCriterion sort_criterion = SCM;
		for(Int4 i=0; i<NUM_LS; i++)
		{
			stable_sort( lattice_result[i].begin(), lattice_result[i].end(), OutputInfo::CmpFunc[(Int4)sort_criterion] );
		}

		// Sets output_flag.
ZLOG_INFO( "Outputting results...\n" );
		OutputInfo outinfo[NUM_LS];
		for(Int4 i=0; i<NUM_LS; i++)
		{
			outinfo[i].setLabel(lattice_result[i], cData);
		}

		// Solution having the top M is output as the best solution.
		printHKLdata(lattice_result, outinfo, sort_criterion, cData, pData,
						cf.putOutputFileName());

ZLOG_INFO( "Input a lattice number in " + cf.putOutputFileName()
		+ "\n(Then, the program outputs an IGOR text file for comparison of calculated peak-positions with the powder diffraction pattern.\n"
		+ "Input \"quit\" to finish the program.) :" );

		string str;
		Int4 cs_index;
		Int4 cs_label;
		vector<LatticeFigureOfMeritToDisplay> selected_lattice_tray;

		do{
			cin >> str;
			if( str == "quit" )
			{
				ZLOG_INFO( "quit.\n" );
				break;
			}

			if( cin.fail() || str.length() < 3
					|| !str2num(str.substr(0,2), cs_index) || !str2num(str.substr(2), cs_label)
					|| cs_index <= 0 || put_number_of_bravais_types() < cs_index 
					|| cs_label <= 0 )
			{
				ZLOG_ERROR( "Wrong lattice number.\n" );
			}
			else
			{
				const Int4 index = outinfo[cs_index-1].putIndex(cs_label);
				
				if( index < 0 )
				{
					ZLOG_ERROR( "Wrong lattice number.\n" );
				}
				else
				{
					const LatticeFigureOfMeritZeroShift& selected_lattice0 = lattice_result[cs_index-1][index].putLatticeFigureOfMerit();
					const size_t n = selected_lattice_tray.size();
					selected_lattice_tray.resize(n+1);

					VecDat3<Double> length_axis, angle_axis;
					selected_lattice0.putReducedLatticeConstantsDegree(cData.putBaseCenteredAxis(), cData.putRhombohedralAxis(), length_axis, angle_axis);

ZLOG_INFO( "Optimizing lattice parameters by linear least squares...\n" );
					bool fitting_succeed = false;
					try{
						// Starts from setting the Bravais-type and unit-cell parameters.
						zerr = selected_lattice_tray[n].setLatticeConstantsDegree(selected_lattice0.enumBravaisType(), cData.putBaseCenteredAxis(), cData.putRhombohedralAxis(), length_axis, angle_axis);
						if( zerr.putErrorType() != ZErrorNoError ) throw zerr;

						selected_lattice_tray[n].putOptimizedLatticeConstantsDegree(cData.putBaseCenteredAxis(), cData.putRhombohedralAxis(), length_axis, angle_axis);
ZLOG_INFO( "Initial unit-cell parameters : "
		+ num2str( length_axis[0] ) + " " + num2str( length_axis[1] ) + " " + num2str( length_axis[2] ) + " "
		+ num2str( angle_axis[0] ) + " " + num2str( angle_axis[1] ) + " " + num2str( angle_axis[2] ) + "\n" );

						zerr = selected_lattice_tray[n].setPeakShiftParamDegree(
															selected_lattice0.putPeakShiftFunctionType(),
															selected_lattice0.putWaveLength(),
															selected_lattice0.putPeakShiftParamDegree(),
															pData);
						if( zerr.putErrorType() != ZErrorNoError ) throw zerr;
						if( selected_lattice0.putPeakShiftFunctionType() != kPeakShiftFunction_Type0 )
						{
ZLOG_INFO( "Wave-length: " + num2str( selected_lattice_tray[n].putLatticeFigureOfMerit().putWaveLength() ) + "\n"
			+ "Initial zero point shift: " + num2str(selected_lattice_tray[n].putPeakShiftParamDegree()[0].value ) + "\n" );
						}

						// Reduce the lattice parameters to be optimized.
						selected_lattice_tray[n].reduceLatticeConstants();

						selected_lattice_tray[n].putOptimizedLatticeConstantsDegree(cData.putBaseCenteredAxis(), cData.putRhombohedralAxis(), length_axis, angle_axis);
ZLOG_INFO( "Reduced unit-cell parameters : "
	+ num2str( length_axis[0] ) + " " + num2str( length_axis[1] ) + " " + num2str( length_axis[2] ) + " "
	+ num2str( angle_axis[0] ) + " " + num2str( angle_axis[1] ) + " " + num2str( angle_axis[2] ) + "\n" );


						// Set figure of merits.
						selected_lattice_tray[n].setFigureOfMerit(cData.putNumberOfReflectionsForFigureOfMerit());

						// Index peaks and set flags.
						selected_lattice_tray[n].resetMillerIndicesInRange(cData.putNumberOfReflectionsForFigureOfMerit());
						selected_lattice_tray[n].resetMillerIndicesToFit();

						// Set Miller indices for optimization.
						zerr = selected_lattice_tray[n].setMillerIndicesToFit(selected_lattice_tray[n].putMillerIndicesToFit());
						if( zerr.putErrorType() != ZErrorNoError ) throw zerr;
						// Set Use/No Use IDs.
						ZErrorMessage zerr = selected_lattice_tray[n].setFittingIDs(selected_lattice_tray[n].putFittingIDs());
						if( zerr.putErrorType() != ZErrorNoError ) throw zerr;

						if( selected_lattice0.putPeakShiftFunctionType() == kPeakShiftFunction_Type0 )
						{
							vector<etype_ID> peak_shift_fitflag(0, _ZRietveldIDVary);
							fitting_succeed = selected_lattice_tray[n].fitLatticeParameter(pData, peak_shift_fitflag, 30, 1.0e-3);	// All arguments are not used in this case.
						}
						else
						{
							vector<etype_ID> peak_shift_fitflag(1, _ZRietveldIDVary);
							fitting_succeed = selected_lattice_tray[n].fitLatticeParameter(pData, peak_shift_fitflag, 30, 1.0e-3);
						}
					}
					catch(const ZErrorMessage& zerr)
					{
ZLOG_ERROR( zerr.printErrorLog() );
					}

					if( fitting_succeed )
					{
						selected_lattice_tray[n].putOptimizedLatticeConstantsDegree(cData.putBaseCenteredAxis(), cData.putRhombohedralAxis(), length_axis, angle_axis);
ZLOG_INFO( "Optimized unit-cell parameters : "
	+ num2str( length_axis[0] ) + " " + num2str( length_axis[1] ) + " " + num2str( length_axis[2] ) + " "
	+ num2str( angle_axis[0] ) + " " + num2str( angle_axis[1] ) + " " + num2str( angle_axis[2] ) + "\n)" );

						if( selected_lattice0.putPeakShiftFunctionType() != kPeakShiftFunction_Type0 )
						{
ZLOG_INFO( "Optimized zero point shift : " + num2str( selected_lattice_tray[n].putPeakShiftParamDegree()[0].value ) + "\n");
						}
						// Save the optimized solution in the original entry in lattice_result.
						lattice_result[cs_index-1][index].setLatticeFigureOfMerit( selected_lattice_tray[n].putLatticeFigureOfMerit() );

						// Reduce the lattice parameters.
						selected_lattice_tray[n].reduceLatticeConstants();
						selected_lattice_tray[n].putOptimizedLatticeConstantsDegree(cData.putBaseCenteredAxis(), cData.putRhombohedralAxis(), length_axis, angle_axis);
ZLOG_INFO( "Reduced unit-cell parameters : "
	+ num2str( length_axis[0] ) + " " + num2str( length_axis[1] ) + " " + num2str( length_axis[2] ) + " "
	+ num2str( angle_axis[0] ) + " " + num2str( angle_axis[1] ) + " " + num2str( angle_axis[2] ) + "\n\n" );
					}

					const Int4 num = selected_lattice_tray[n].putLatticeFigureOfMerit().checkDominantZone();
					if( num > cData.putNumberOfReflectionsForFigureOfMerit() )
					{
						ZLOG_WARN( "A dominant zone exists.  This is resolved by increasing the parameter <" + cData.putNumberOfReflectionsForFigureOfMeritLabel() + "> to a number more than " + num2str(num-1) + ".\n" );
					}

					ZLOG_INFO( "Checking coincidence of computed lines ...\n" );
					LatticeWithSameQ lws;
					lws.setParam( cData, selected_lattice_tray[n].putLatticeFigureOfMerit().putSellingReducedForm() );
					pair<bool, Double> evaluation_ans;
					lws.execute( evaluation_ans, selected_lattice_tray[n].putLatticeFigureOfMerit().putSellingReducedForm() );

					if( !evaluation_ans.first )
					{
						ZLOG_WARN( "it is necessary to increase \"" + cData.putMaxQTocheckComputedLinesLabel()
									+ "\" or reduce " + cData.putResolutionTocheckComputedLinesLabel()
									+ "\" in order to obtain all the unit-cells with the same computed lines.\n" );
					}

					vector<LatticeFigureOfMerit> lattices_same_q;
					lws.putLatticesWithSameQ(cData, selected_lattice_tray[n].putQDataModifiedWithNewPeakShiftParam(), lattices_same_q);

					if( lattices_same_q.empty() )
					{
						ZLOG_INFO( "The unit-cell does not have the same computed lines as any other unit-cells.\n\n" );
					}
					else
					{
				      	ostringstream oss;
						printSolutions(lattices_same_q, cData, &oss);
						ZLOG_WARN( "The unit-cell have very similar computed lines as the following (the solution might not determined uniquely from peak positions):\n" + oss.str() );
						if( lattices_same_q.size() > 15 )
						{
							ZLOG_WARN( "The peak positions of only the first 15 unit-cells are displayed. Output the IGOR file to check all the unit cells." );
						}
					}

					string fname00, fname0;
			    	removeFileExtension(cf.putOutputFileName(), fname00);
			    	removeFileExtension(fname00, fname0);

					const string output_igor_file_name = fname0 + "_lattice("
							+ put_bravais_type_name(selected_lattice_tray[n].putLatticeFigureOfMerit().enumBravaisType(), cData.putBaseCenteredAxis()) + ";"
							+ selected_lattice_tray[n].putLatticeFigureOfMerit().printOptimizedLatticeConstants(cData.putBaseCenteredAxis(), cData.putRhombohedralAxis(), 3) + ";"
							+ num2str<Double>(selected_lattice_tray[n].putLatticeFigureOfMerit().putFiguresOfMerit().putFigureOfMeritWolff(), 3) + ").histogramIgor";

					printPeakPosition(cData, pData, selected_lattice_tray[n], lattices_same_q, output_igor_file_name);

			  		ZLOG_INFO( "The program has output an IGOR text file : " + output_igor_file_name + "\n\n" );

			        if( false )
			        {
						const Int4 num_types = putNumberOfTypesOfSystematicAbsences(selected_lattice_tray[n].putLatticeFigureOfMerit().putBravaisType());
						const Double MIN_FOM = floor(selected_lattice_tray[n].putLatticeFigureOfMerit().putFiguresOfMerit().putFigureOfMeritWolff()*100.0)*0.01;

						vector<LatticeFigureOfMeritToDisplay> lattices_with_reflection_condition;
				        for(Int4 i=0; i<num_types; i++)
				        {
				            // Copy constructor.
				        	LatticeFigureOfMeritToDisplay lattice_wrc(selected_lattice_tray[n]);

				        	lattice_wrc.setTypeOfSystematicAbsences(i);
				           	lattice_wrc.setDeWolffFigureOfMerit(cData.putNumberOfReflectionsForFigureOfMerit());
				    		if( lattice_wrc.putLatticeFigureOfMerit().putFiguresOfMerit().putFigureOfMeritWolff() < MIN_FOM ) continue;

				    		lattice_wrc.resetMillerIndicesInRange(cData.putNumberOfReflectionsForFigureOfMerit());
				    		lattice_wrc.resetMillerIndicesToFit();

				            lattices_with_reflection_condition.push_back(lattice_wrc);
				        }
				        sort( lattices_with_reflection_condition.begin(), lattices_with_reflection_condition.end(), cmpDeWolff );

		//		        print_lattice_array(cData, pData, lattices_with_reflection_condition, fname0 + ".out.xml");

						const string output_igor_file_name = fname0 + "_lattice_ref_cond("
								+ put_bravais_type_name(selected_lattice_tray[n].putLatticeFigureOfMerit().enumBravaisType(), cData.putBaseCenteredAxis()) + ";"
								+ selected_lattice_tray[n].putLatticeFigureOfMerit().printOptimizedLatticeConstants(cData.putBaseCenteredAxis(), cData.putRhombohedralAxis(), 3) + ";"
								+ num2str<Double>(selected_lattice_tray[n].putLatticeFigureOfMerit().putFiguresOfMerit().putFigureOfMeritWolff(), 3) + ").histogramIgor";

						printPeakPosition(cData, pData, MIN_FOM,
								lattices_with_reflection_condition, output_igor_file_name);
			        }
				}
			}

			ZLOG_INFO( "Input a lattice number in " + cf.putOutputFileName() + " :" );

  		} while( true );

		// Selected lattice is output.
		if( !selected_lattice_tray.empty() )
		{
			string fname00, fname0;
	    	removeFileExtension(cf.putOutputFileName(), fname00);
	    	removeFileExtension(fname00, fname0);

	    	printHKLdata(selected_lattice_tray, cData, pData, fname0 +".index2.xml");

			printPeakPosition(cData, pData, selected_lattice_tray, fname0 +"_lattices.histogramIgor");
		}
	}
	catch(bad_alloc& ball){
		ZErrorMessage zerr = nerror(ball, __FILE__, __LINE__, __FUNCTION__);
		ZLOG_ERROR( zerr.printErrorLog() );
		return 0;
	}
	catch(out_of_range&)
	{
		ZErrorMessage zerr("out_of_range exception has occurred", __FILE__, __LINE__, __FUNCTION__);
		ZLOG_ERROR( zerr.printErrorLog() );
		return 0;
	}
	catch(const ZErrorMessage& etype)
	{
		ZLOG_ERROR( etype.printErrorLog() );
		return 0;
	}

    return 1;
}

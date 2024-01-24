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
#include <iostream>
#include "p_out_indexing.hh"
#include "ControlParam.hh"
#include "PeakPosData.hh"
#include "lattice_symmetry/VCLatticeFigureOfMeritToCheckSymmetry.hh"
#include "lattice_symmetry/LatticeFigureOfMeritToDisplay.hh"
#include "utility_func/covar_matrix.hh"
#include "utility_func/zmath.hh"
#include "utility_lattice_reduction/matrix_NbyN.hh"
#include "utility_data_structure/Bud.hh"
#include "utility_data_structure/Node3.hh"
#include "utility_data_structure/qdata.hh"
#include "utility_data_structure/SymMatWCovar.hh"
#include "utility_data_structure/VecDat3.hh"
#include "lattice_symmetry/OutputInfo.hh"


void printQdata(const vector<QData>& qdata,
		ostream* os)
{
    os->setf(ios::right);
    os->setf(ios::uppercase);
    os->setf(ios::showpoint);
    os->setf(ios::scientific);

    *os << "IGOR" << endl;
    *os << "WAVES/O index, Q, sqrt_root_of_variance_of_Q" << endl;
    *os << "BEGIN" << endl;

    const Int4 isize = qdata.size();
    
    for (Int4 k=0; k<isize; k++)
    {
       	os->width(5);
       	*os << k+1;

       	os->width(15);
       	*os << qdata[k].q;

       	os->width(15);
       	*os << sqrt(qdata[k].q_var);

       	*os << endl;
    }
    *os << "END" << endl;
}



void print_bud_data(const vector<Bud>& bud_basket, 
const Vec_BOOL& root_flag, 
const vector< vector<Int4> >& leftbr_tray,
const vector< vector<Int4> >& rightbr_tray,
const vector< vector<Int4> >& root_leftbr,
const vector< vector<Int4> >& root_rightbr,
const string& fname)
{
    ofstream ofs(fname.c_str());
    ostream *os;
    os = &ofs;
    
    *os << "IGOR" << endl;
    *os << "WAVES/O ";
    
    *os << "Number, K1, K2, K3, K4, (Q1+Q2)*2-Q3-Q4, Err[(Q1+Q2)*2-Q3-Q4]";
    *os << "root?, super_base?, left_branch, right_branch, left_root, right_root, tree_index" << endl;
    
    *os << "BEGIN" << endl;
    os->setf(ios::right);
    os->setf(ios::uppercase);
    os->setf(ios::showpoint);

    Int4 k=0;
    for(vector<Bud>::const_iterator it=bud_basket.begin(); it!=bud_basket.end(); it++, k++)
    {
       	os->unsetf(ios::scientific);
       	os->width(5);
       	*os << k+1;

       	os->width(5);
       	if( it->iK1() < 0 ) *os << -1;
       	else *os << it->iK1() + 1;

       	os->width(5);
       	if( it->iK2() < 0 ) *os << -1;
       	else *os << it->iK2() + 1;

       	os->width(5);
       	if( it->iK3() < 0 ) *os << -1;
       	else *os << it->iK3() + 1;

       	os->width(5);
       	if( it->iK4() < 0 ) *os << -1;
       	else *os << it->iK4() + 1;

       	const VCData diff = (it->Q1()+it->Q2())*2-(it->Q3()+it->Q4());

       	os->setf(ios::scientific);
       	os->width(15);
       	*os << diff.Value();

       	os->width(15);
       	*os << sqrt( diff.Variance() );

       	os->precision(6);
       	os->width(5);
    	if( root_flag[k] ) *os << 1;
    	else *os << 0;
    		
       	os->width(5);
    	if( it->IsSuperBasisObtuse() ) *os << 1;
    	else *os << 0;
    		
       	os->width(10);
       	if( leftbr_tray[k].empty() ) *os << "Emp";
    	else
    	{
    		*os << leftbr_tray[k][0] + 1;
        	for(UInt4 l=1; l<leftbr_tray[k].size(); l++) *os << ":" << leftbr_tray[k][l] + 1;
    	}

       	os->width(10);
    	if( rightbr_tray[k].empty() ) *os << "Emp";
    	else
    	{
    		*os << rightbr_tray[k][0] + 1;
        	for(UInt4 l=1; l<rightbr_tray[k].size(); l++) *os << ":" << rightbr_tray[k][l] + 1;
    	}

       	os->width(10);
    	if( root_leftbr[k].empty() ) *os << "Emp";
    	else
    	{
    		*os << root_leftbr[k][0] + 1;
        	for(UInt4 l=1; l<root_leftbr[k].size(); l++) *os << ":" << root_leftbr[k][l] + 1;
    	}

       	os->width(10);
    	if( root_rightbr[k].empty() ) *os << "Emp";
    	else
    	{
    		*os << root_rightbr[k][0] + 1;
        	for(UInt4 l=1; l<root_rightbr[k].size(); l++) *os << ":" << root_rightbr[k][l] + 1;
    	}
    	
    	*os << endl;
    }
    *os << "END" << endl;

    ofs.close();
}




void print_bud_data(const vector<Bud>& bud_basket, const string& fname)
{
    ofstream ofs(fname.c_str());
    ostream *os;
    os = &ofs;

    *os << "IGOR" << endl;
    *os << "WAVES/O ";

    *os << "Number, K1, K2, K3, K4, (Q1+Q2)*2-Q3-Q4, Err[(Q1+Q2)*2-Q3-Q4]";
    *os << "root?, super_base?, left_branch, right_branch, left_root, right_root, tree_index" << endl;

    *os << "BEGIN" << endl;
    os->setf(ios::right);
    os->setf(ios::uppercase);
    os->setf(ios::showpoint);

    Int4 k=0;
    for(vector<Bud>::const_iterator it=bud_basket.begin(); it!=bud_basket.end(); it++, k++)
    {
       	os->unsetf(ios::scientific);
       	os->width(5);
       	*os << k+1;

       	os->width(5);
       	if( it->iK1() < 0 ) *os << -1;
       	else *os << it->iK1() + 1;

       	os->width(5);
       	if( it->iK2() < 0 ) *os << -1;
       	else *os << it->iK2() + 1;

       	os->width(5);
       	if( it->iK3() < 0 ) *os << -1;
       	else *os << it->iK3() + 1;

       	os->width(5);
       	if( it->iK4() < 0 ) *os << -1;
       	else *os << it->iK4() + 1;

       	const VCData diff = (it->Q1()+it->Q2())*2-(it->Q3()+it->Q4());

       	os->setf(ios::scientific);
       	os->width(15);
       	*os << diff.Value();

       	os->width(15);
       	*os << sqrt( diff.Variance() );

       	os->width(5);
    	if( it->IsSuperBasisObtuse() ) *os << 1;
    	else *os << 0;

    	*os << endl;
    }
    *os << "END" << endl;

    ofs.close();
}




void print_node_data(const vector<Node3>& node_basket,
					const string& fname)
{
	ofstream ofs(fname.c_str());
    ostream *os;
    os = &ofs;

	ofs << "** Node, K1, K2, K3, K1-K2, K1-K3, K2-K3, Q1, Q2, Q3, Q12, Q13, Q23, UnitcellVolume" << endl;

   	os->setf(ios::scientific);
   	os->precision(4);

   	Int4 iK1, iK2, iK3;
   	Int4 iK12, iK13, iK23;
	Double detS;
   	for(UInt4 l=0; l<node_basket.size(); l++)
	{
       	os->width(5);
       	*os << l + 1;

       	const Node3& nodex = node_basket[l];
       	iK1 = nodex.iK(0);
       	iK2 = nodex.iK(1);
       	iK3 = nodex.iK(2);

       	iK12 = nodex.iSumK(0,1);
       	iK13 = nodex.iSumK(0,2);
       	iK23 = nodex.iSumK(1,2);

       	os->width(8);
      	if( iK1 >= 0 ) *os << iK1 + 1;
      	else *os << -1;

       	os->width(5);
      	if( iK2 >= 0 ) *os << iK2 + 1;
      	else *os << -1;

       	os->width(5);
      	if( iK3 >= 0 ) *os << iK3 + 1;
      	else *os << -1;

      	os->width(8);
      	if( iK12 >= 0 ) *os << iK12 + 1;
      	else *os << -1;

       	os->width(5);
      	if( iK13 >= 0 ) *os << iK13 + 1;
      	else *os << -1;

       	os->width(5);
      	if( iK23 >= 0 ) *os << iK23 + 1;
      	else *os << -1;


      	os->width(15);
      	*os << nodex.K(0).Value();

       	os->width(15);
      	*os << nodex.K(1).Value();

       	os->width(15);
      	*os << nodex.K(2).Value();

      	os->width(15);
      	*os << nodex.SumK(0,1).Value();

       	os->width(15);
      	*os << nodex.SumK(0,2).Value();

       	os->width(15);
      	*os << nodex.SumK(1,2).Value();

		detS = nodex.putDeterminantS();

      	os->width(20);
       	*os << 1.0/sqrt(detS);

       	*os << endl;
	}
    *os << endl;
}




void printHKLdata(const vector<VCLatticeFigureOfMeritToCheckSymmetry> lattice_result[],
		const OutputInfo outinfo[],
		const eSortCriterion& sort_criterion,
		const ControlParam& cdata,
		const PeakPosData& pdata,
		const string& fname)
{
	static const Int4 NUM_LS = put_number_of_bravais_types();

	ofstream ofs(fname.c_str());
    ostream *os;
    os = &ofs;
    
   	os->setf(ios::scientific);
   	os->precision(4);
   	
	pair< vector<Int4>::const_iterator, vector<Int4>::const_iterator> it_pair;
	SymMatWCovar dbl_S(3);
	VecDat3<Double> length_axis, angle_axis; 
	SymMat<VCData> S_red2(3);
	set< SymMat<VCData> > S_tray;

	Int4 label_start=0;
	os->width(label_start);
	*os << "" << "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\n";

	os->width(label_start);
	*os << "" << "<ZCodeParameters>\n";
	label_start++;

	os->width(label_start);
	*os << "" << "<ConographOutput>\n";
	label_start++;
	
	os->width(label_start);
	*os << "" << "<!-- Information on the best " + putLabel(sort_criterion) + " solution for each Bravais type.\n";
	os->width(label_start+6);
	*os << "" << "TNB : number of solutions with the Bravais type,\n";
	os->width(label_start+6);
	*os << "" << putLabel(SCM) << " : de Wolff figure of merit,\n";
	os->width(label_start+6);
	*os << "" << putLabel(SCWuM) << " : Wu figure of merit,\n";
//	os->width(label_start+6);
//	*os << "" << putLabel(SCNormM) << " : normalized de Wolff figure of merit,\n";
	os->width(label_start+6);
	*os << "" << putLabel(SCRevM) << " : reversed de Wolff figure of merit,\n";
	os->width(label_start+6);
	*os << "" << putLabel(SCSymM) << " : symmetric de Wolff figure of merit,\n";
	os->width(label_start+6);
	*os << "" << putLabel(SCNN) << " : number of lattices in the neighborhood,\n";
	os->width(label_start+6);
	*os << "" << "VOL : unit-cell volume.\n";
	os->width(label_start);
	*os << "" << "Bravais Lattice : TNB";
	for(Int4 i=0; i<putNumberOfSortCriterion(); i++)
	{
		*os << ", " << putLabel(eSortCriterion(i));
	}
	*os << ", VOL\n";

	for(Int4 k=NUM_LS-1; k>=0; k--)
	{
		os->width(17);
		*os << put_bravais_type_name(eBravaisType(k), cdata.putBaseCenteredAxis()) << " : ";
		os->width(7);
		*os << outinfo[k].putCandidateNumToOutput();

		const Int4 index = outinfo[k].putIndex( outinfo[k].putLatticeLabelSelectedFromListToOutput(sort_criterion) );
		if( index >= 0 )
		{
			const LatticeFigureOfMerit& lat_fom = lattice_result[k][index].putLatticeFigureOfMerit();
			const LatticeFigureOfMerit::SetOfFigureOfMerit& set_fom = lat_fom.putFiguresOfMerit();

			os->width(15);
			*os << set_fom.putFigureOfMeritWolff();

			os->width(15);
			*os << set_fom.putFigureOfMeritWu();

//			os->width(15);
//			*os << set_fom.putNormalizedFigureOfMeritWolff();

			os->width(15);
			*os << set_fom.putReversedFigureOfMerit();

	   	   	os->width(15);
			*os << set_fom.putSymmetricFigureOfMerit();

			os->width(5);
			*os << lattice_result[k][index].putNumberOfLatticesInNeighborhood();

			os->width(15);
			*os << lat_fom.putLatticeVolume();

//			os->width(15);
//			*os << set_fom.putFigureOfMeritWolff_Original();
//
//			os->width(15);
//			*os << set_fom.putFigureOfMeritWu_Original();
		}
		*os << endl;
	}
	os->width(label_start);
	*os << "" << "-->\n\n";
	os->width(label_start);
	*os << "" << "<!-- Labels of the solution with the best figure of merit.\n";
	for(Int4 k=NUM_LS-1; k>=0; k--)
	{
		os->width(18);
		*os << put_bravais_type_name(eBravaisType(k), cdata.putBaseCenteredAxis()) << ", Best Score : Lattice constants, label.\n";
//		if( lattice_result[k].empty() )
//		{
//			*os << "\n";
//			continue;
//		}
		for(Int4 i=0; i<putNumberOfSortCriterion(); i++)
		{
			os->width(15);
			*os << putLabel(eSortCriterion(i)) << " = ";

			const Int4 index = outinfo[k].putIndex( outinfo[k].putLatticeLabelSelectedAmongAll(eSortCriterion(i)) );
			os->width(12);
			if( index >= 0 )
			{
				const LatticeFigureOfMerit& lat_fom = lattice_result[k][index].putLatticeFigureOfMerit();
				const LatticeFigureOfMerit::SetOfFigureOfMerit& set_fom = lat_fom.putFiguresOfMerit();

				if( eSortCriterion(i) == SCM )
				{
					*os << set_fom.putFigureOfMeritWolff();
				}
				else if( eSortCriterion(i) == SCWuM )
				{
					*os << set_fom.putFigureOfMeritWu();
				}
//				if( eSortCriterion(i) == SCNormM )
//				{
//					*os << set_fom.putNormalizedFigureOfMeritWolff();
//				}
				else if( eSortCriterion(i) == SCRevM )
				{
					*os << set_fom.putReversedFigureOfMerit();
				}
				else if( eSortCriterion(i) == SCSymM )
				{
					*os << set_fom.putSymmetricFigureOfMerit();
				}
				else if( eSortCriterion(i) == SCNN )
				{
					*os << lattice_result[k][index].putNumberOfLatticesInNeighborhood();
				}
//				else if( eSortCriterion(i) == SCM_org )
//				{
//					*os << set_fom.putFigureOfMeritWolff_Original();
//				}
//				else if( eSortCriterion(i) == SCWuM_org )
//				{
//					*os << set_fom.putFigureOfMeritWu_Original();
//				}

				*os << " : ";
				lat_fom.putReducedLatticeConstantsDegree(cdata.putBaseCenteredAxis(), cdata.putRhombohedralAxis(), length_axis, angle_axis);
			  	os->width(14);
			  	*os << length_axis[0];
			  	os->width(14);
			   	*os << length_axis[1];
			  	os->width(14);
			   	*os << length_axis[2];
			 	os->width(14);
			   	*os << angle_axis[0];
			  	os->width(14);
			   	*os << angle_axis[1];
			  	os->width(14);
			   	*os << angle_axis[2];
			  	os->width(10);
			  	*os << lattice_result[k][index].putStringLabel();
			  	*os << endl;
			}
			else
			{
				*os << "- -" << " : - -\n";
			}
		}
		*os << endl;
	}
	os->width(label_start);
	*os << "" << "-->\n\n";

	pair<eBravaisType, lattice_label> selected_lattice_info = put_selected_lattice_from_all(lattice_result, outinfo, sort_criterion);
	const Int4 selected_lattice_index = outinfo[selected_lattice_info.first].putIndex(selected_lattice_info.second);

	if( selected_lattice_index >= 0 )
	{
		const VCLatticeFigureOfMeritToCheckSymmetry& selected_lattice = lattice_result[selected_lattice_info.first][selected_lattice_index];

		os->width(label_start);
		*os << "" << "<!-- Information on the selected candidates.-->\n";

		os->width(label_start);
	  	*os << "" << "<SelectedLatticeCandidate number=\"" << selected_lattice.putStringLabel() << "\">\n";
	  	label_start++;

	  	selected_lattice.printLatticeInformation(lattice_result, outinfo,
	  										cdata.putBaseCenteredAxis(), cdata.putRhombohedralAxis(),
	  										cdata.putResolution(), label_start, os);

	  	label_start--;
		os->width(label_start);
	  	*os << "" << "</SelectedLatticeCandidate>\n\n";
	}
	
	for(Int4 k=NUM_LS-1; k>=0; k--)
	{
		*os << "\n";
		os->width(label_start);
		*os << "" << "<!-- Candidates for " << put_bravais_type_name(eBravaisType(k), cdata.putBaseCenteredAxis()) << " -->\n\n";
		
		const Int4 num_topo = lattice_result[k].size();

		for(Int4 n=0; n<num_topo; n++)
	   	{
			const VCLatticeFigureOfMeritToCheckSymmetry& ans = lattice_result[k][n];
			if( !outinfo[k].IsOutput(ans.putLatticeLabel()) ) continue;

			os->width(label_start);
	      	*os << "" << "<LatticeCandidate number=\"" << ans.putStringLabel() << "\">\n";
	      	label_start++;
	      	ans.printLatticeInformation(lattice_result, outinfo, cdata.putBaseCenteredAxis(), cdata.putRhombohedralAxis(),
	      										0.06, label_start, os);
	      	
	      	label_start--;
			os->width(label_start);
	      	*os << "" << "</LatticeCandidate>\n\n";
	   	}
   	   	*os << endl;
	}
	
	label_start--;
	os->width(label_start);
	*os << "" << "</ConographOutput>\n";

	label_start--;
	os->width(label_start);
	*os << "" << "</ZCodeParameters>\n";
}



void printHKLdata(const vector<LatticeFigureOfMeritToDisplay> lattice_result,
		const ControlParam& cdata,
		const PeakPosData& pdata,
		const string& fname)
{
	ofstream ofs(fname.c_str());
    ostream *os;
    os = &ofs;

   	os->setf(ios::scientific);
   	os->precision(4);

	pair< vector<Int4>::const_iterator, vector<Int4>::const_iterator> it_pair;
	SymMatWCovar dbl_S(3);
	VecDat3<Double> length_axis, angle_axis;
	SymMat<VCData> S_red2(3);
	set< SymMat<VCData> > S_tray;

	Int4 label_start=0;
	os->width(label_start);
	*os << "" << "<ZCodeParameters>\n";
	label_start++;

	os->width(label_start);
	*os << "" << "<ConographOutput>\n";
	label_start++;

	const Int4 num_topo = lattice_result.size();

	for(Int4 n=0; n<num_topo; n++)
   	{
		const LatticeFigureOfMeritZeroShift& ans = lattice_result[n].putLatticeFigureOfMerit();

		os->width(label_start);
      	*os << "" << "<LatticeCandidate>\n";
      	label_start++;
      	ans.printLatticeInformation(cdata.putBaseCenteredAxis(), cdata.putRhombohedralAxis(),
      								cdata.putResolution(), label_start, os);
      	lattice_result[n].printIndexingResult(cdata, pdata, label_start, os);

      	label_start--;
		os->width(label_start);
      	*os << "" << "</LatticeCandidate>\n\n";
   	}
   	*os << endl;

	label_start--;
	os->width(label_start);
	*os << "" << "</ConographOutput>\n";

	label_start--;
	os->width(label_start);
	*os << "" << "</ZCodeParameters>\n";
}




void printPeakPosition(
		const ControlParam& cdata,
		const PeakPosData& pdata,
		const LatticeFigureOfMeritToDisplay& latfit,
		const vector<LatticeFigureOfMerit>& lattices_same_q,
		const string& fname)
{
	ofstream ofs(fname.c_str());
    ostream * const os = &ofs;

    *os << "IGOR\n";
    *os << "WAVES/O ";

    pdata.printData(os);

	*os << "WAVES/O dphase_0, xphase_0, yphase_0, h_0, k_0, l_0\n";

    *os << "BEGIN\n";

	const vector<HKL_Q>& cal_hkl_tray = latfit.putCalMillerIndices();
	Vec_DP cal_pos_tray;
	latfit.putCalculatedPeakPosInRange(cdata, cal_pos_tray);

	const Int4 peak_num = cal_hkl_tray.size();

	for (Int4 i=0; i<peak_num; i++)
	{
        os->width(15);
        *os << 1.0 / sqrt( cal_hkl_tray[i].Q() );

        os->width(15);
        if( cal_pos_tray[i] < 0.0 )
        {
            *os << "NAN";
        }
        else
        {
        	*os << cal_pos_tray[i];
        }
        os->width(15);
        *os << 0.0;

        const VecDat3<Int4>& hkl = cal_hkl_tray[i].HKL();

        os->width(5);
        *os << hkl[0];
        os->width(5);
        *os << hkl[1];
        os->width(5);
        *os << hkl[2];

        *os << endl;
	}
	*os << "END\n\n";

	LatticeFigureOfMeritToDisplay latfit2;
	latfit2.setPeakShiftParamDegree(latfit.putPeakShiftFunctionType(), latfit.putWaveLength(),
										latfit.putPeakShiftParamDegree(), pdata);

	const Int4 isize = lattices_same_q.size();
    for (Int4 j=0, j2=1; j<isize; j++, j2++)
    {
		VecDat3<Double> length_axis, angle_axis;
		lattices_same_q[j].putReducedLatticeConstantsDegree(cdata.putBaseCenteredAxis(), cdata.putRhombohedralAxis(), length_axis, angle_axis);

		latfit2.setLatticeConstantsDegree(lattices_same_q[j].enumBravaisType(), cdata.putBaseCenteredAxis(), cdata.putRhombohedralAxis(), length_axis, angle_axis);
		latfit2.resetMillerIndicesInRange(cdata.putNumberOfReflectionsForFigureOfMerit());


    	*os << "WAVES/O dphase_" << j2 << ", xphase_" << j2 << ", yphase_" << j2 << ", ";
    	*os << "h_" << j2 << ", ";
    	*os << "k_" << j2 << ", ";
    	*os << "l_" << j2;
    	*os << endl;

	    *os << "BEGIN\n";

		const vector<HKL_Q>& cal_hkl_tray = latfit2.putCalMillerIndices();
		Vec_DP cal_pos_tray;
		latfit2.putCalculatedPeakPosInRange(cdata, cal_pos_tray);

    	const Int4 peak_num = cal_hkl_tray.size();

    	for (Int4 i=0; i<peak_num; i++)
    	{
            os->width(15);
            *os << 1.0 / sqrt( cal_hkl_tray[i].Q() );

            os->width(15);
            if( cal_pos_tray[i] < 0.0 )
            {
                *os << "NAN";
            }
            else
            {
            	*os << cal_pos_tray[i];
            }
            os->width(15);
            *os << 0.0;

            const VecDat3<Int4>& hkl = cal_hkl_tray[i].HKL();

            os->width(5);
            *os << hkl[0];
            os->width(5);
            *os << hkl[1];
            os->width(5);
            *os << hkl[2];

            *os << endl;
    	}
    	*os << "END\n\n";
    }

    VecDat3<Double> length_axis, angle_axis;
    const string str_num_ref = num2str( cdata.putNumberOfReflectionsForFigureOfMerit() );


    *os << "X Display " << pdata.putYIntColumnTitle() << " vs " << pdata.putXColumnTitle() << endl;
	*os << "X ModifyGraph mirror(left)=2\n";
	*os << "X ModifyGraph mirror(bottom)=2\n";
	*os << "X ModifyGraph rgb(" << pdata.putYIntColumnTitle() << ")=(0,15872,65280)\n";

    *os << "X AppendToGraph height vs peakpos\n";
    *os << "X ModifyGraph mode(height)=3,marker(height)=17\n";
    *os << "X ModifyGraph rgb(height)=(0,65280,65280)\n";

    Double offset = 0.0;
    const Double offset_gap = pdata.putMaxPeakHeightOfFirst20() * (-0.05);
    for (Int4 j=0; j<isize+1; j++, offset+=offset_gap)
	{
		*os << "X AppendToGraph yphase_" << j << " vs xphase_" << j << endl;
		*os << "X ModifyGraph offset(yphase_" << j << ")={0," << offset << "},mode(yphase_" << j << ")=3,marker(yphase_" << j;
		*os << ")=10,msize(yphase_" << j << ")=3,mrkThick(yphase_" << j << ")=0.6,rgb(yphase_" << j << ")=(3,52428,1)" << endl;
    	if( j <= 0 ) offset+=offset_gap;
	}

    *os << "X Legend/C/N=text0/J/A=MC \"";
	*os << "\\s(yphase_" << 0 << ") "
		+ put_bravais_type_name(latfit.putLatticeFigureOfMerit().enumBravaisType(), cdata.putBaseCenteredAxis()) + " "
		+ latfit.putLatticeFigureOfMerit().printOptimizedLatticeConstants(cdata.putBaseCenteredAxis(), cdata.putRhombohedralAxis(), 3);
	for (Int4 j=0, j2=1; j<isize; j++, j2++)
	{
		*os << "\\r\\s(yphase_" << j2 << ") "
				+ put_bravais_type_name(lattices_same_q[j].enumBravaisType(), cdata.putBaseCenteredAxis()) + " "
				+ lattices_same_q[j].printOptimizedLatticeConstants(cdata.putBaseCenteredAxis(), cdata.putRhombohedralAxis(), 3);
	}
	*os << "\"\nX Legend/C/N=text0/J/A=RT/X=0.00/Y=0.00\n";

    ofs.close();
}



void printPeakPosition(
		const ControlParam& cdata,
		const PeakPosData& pdata,
		const vector<LatticeFigureOfMeritToDisplay>& latfit_tray,
		const string& fname)
{
	ofstream ofs(fname.c_str());
    ostream * const os = &ofs;

    *os << "IGOR\n";
    *os << "WAVES/O ";

    pdata.printData(os);

    if( latfit_tray.empty() ) return;

    const Int4 isize = latfit_tray.size();
    for (Int4 j=0, j1=1; j<isize; j++, j1++)
    {
    	const LatticeFigureOfMeritToDisplay& latfit = latfit_tray[j];

    	*os << "WAVES/O dphase_" << j1 << ", xphase_" << j1 << ", yphase_" << j1 << ", ";
    	*os << "h_" << j1 << ", ";
    	*os << "k_" << j1 << ", ";
    	*os << "l_" << j1;
    	*os << endl;

	    *os << "BEGIN\n";

		const vector<HKL_Q>& cal_hkl_tray = latfit.putCalMillerIndices();
		Vec_DP cal_pos_tray;
		latfit.putCalculatedPeakPosInRange(cdata, cal_pos_tray);

    	const Int4 peak_num = cal_hkl_tray.size();

    	for (Int4 i=0; i<peak_num; i++)
    	{
            os->width(15);
            *os << 1.0 / sqrt( cal_hkl_tray[i].Q() );

            os->width(15);
            if( cal_pos_tray[i] < 0.0 )
            {
                *os << "NAN";
            }
            else
            {
            	*os << cal_pos_tray[i];
            }
            os->width(15);
            *os << 0.0;

            const VecDat3<Int4>& hkl = cal_hkl_tray[i].HKL();

            os->width(5);
            *os << hkl[0];
            os->width(5);
            *os << hkl[1];
            os->width(5);
            *os << hkl[2];

            *os << endl;
    	}
    	*os << "END\n\n";
    }

    VecDat3<Double> length_axis, angle_axis;
    const string str_num_ref = num2str( cdata.putNumberOfReflectionsForFigureOfMerit() );


    *os << "X Display " << pdata.putYIntColumnTitle() << " vs " << pdata.putXColumnTitle() << endl;
	*os << "X ModifyGraph mirror(left)=2\n";
	*os << "X ModifyGraph mirror(bottom)=2\n";
	*os << "X ModifyGraph rgb(" << pdata.putYIntColumnTitle() << ")=(0,15872,65280)\n";

    *os << "X AppendToGraph height vs peakpos\n";
    *os << "X ModifyGraph mode(height)=3,marker(height)=17\n";
    *os << "X ModifyGraph rgb(height)=(0,65280,65280)\n";

    Double offset = 0.0;
    const Double offset_gap = pdata.putMaxPeakHeightOfFirst20() * (-0.05);
    for (Int4 j=0, j1=1; j<isize; j++, j1++, offset+=offset_gap)
	{
		*os << "X AppendToGraph yphase_" << j1 << " vs xphase_" << j1 << endl;
		*os << "X ModifyGraph offset(yphase_" << j1 << ")={0," << offset << "},mode(yphase_" << j1 << ")=3,marker(yphase_" << j1;
		*os << ")=10,msize(yphase_" << j1 << ")=3,mrkThick(yphase_" << j1 << ")=0.6,rgb(yphase_" << j1 << ")=(3,52428,1)" << endl;
	}

    *os << "X Legend/C/N=text0/J/A=MC \"";
	*os << "\\s(yphase_" << 1 << ") "
		+ put_bravais_type_name(latfit_tray[0].putLatticeFigureOfMerit().enumBravaisType(), cdata.putBaseCenteredAxis()) + " "
		+ latfit_tray[0].putLatticeFigureOfMerit().printOptimizedLatticeConstants(cdata.putBaseCenteredAxis(), cdata.putRhombohedralAxis(), 3);
	for (Int4 j=1, j1=2; j<isize; j++, j1++)
	{
		*os << "\\r\\s(yphase_" << j1 << ") "
				+ put_bravais_type_name(latfit_tray[j].putLatticeFigureOfMerit().enumBravaisType(), cdata.putBaseCenteredAxis()) + " "
				+ latfit_tray[j].putLatticeFigureOfMerit().printOptimizedLatticeConstants(cdata.putBaseCenteredAxis(), cdata.putRhombohedralAxis(), 3);
	}
	*os << "\"\nX Legend/C/N=text0/J/A=RT/X=0.00/Y=0.00\n";

    ofs.close();
}

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
#include <fstream>
#include "p_out_space_group_dtm.hh"
#include "../ControlParam.hh"
#include "../PeakPosData.hh"
#include "../lattice_symmetry/LatticeFigureOfMeritToDisplay.hh"


void print_lattice_array(const ControlParam& cdata,
		const PeakPosData& pdata,
		const vector<LatticeFigureOfMeritToDisplay>& arg, const string& fname)
{
	ofstream ofs(fname.c_str());

    ofs.setf(ios::right);
    ofs.setf(ios::scientific);
    ofs.precision(6);

    Int4 label_start=0;

    ofs.width(label_start);
    ofs << "" << "<ZCodeParameters>\n";
    label_start++;

    ofs.width(label_start);
    ofs << "" << "<ConographOutput>\n";
    label_start++;

    ofs.width(label_start);
    ofs << "" << "<TypeOfReflectionConditions>\n";
    //label_start++;

	for (size_t j=0; j<arg.size(); j++)
	{
        ofs.width(label_start+1);
        ofs << "" << "<Candidates>\n";
        //label_start++;

        ofs.width(label_start+2);
        ofs << "" << "<SpaceGroups> "+arg[j].putStringTypeOfSystematicAbsences();
        ofs << "" << " </SpaceGroups>\n";
        //label_start++;

        ofs.width(label_start+2);
        ofs << "" << "<ReflectionConditions> "+arg[j].putStringReflectionConditions();
        ofs << "" << " </ReflectionConditions>\n";
        //label_start++;


        ofs.width(label_start+2); ofs << "";
        ofs << "<FigureOfMeritWolff name=\""+arg[j].putLatticeFigureOfMerit().putFiguresOfMerit().putLabel_FigureOfMeritWolff() <<  "\"> ";
        ofs << arg[j].putLatticeFigureOfMerit().putFiguresOfMerit().putFigureOfMeritWolff() << " </FigureOfMeritWolff>\n";

      	arg[j].printIndexingResult(cdata, pdata, label_start+2, &ofs);

      	ofs.width(label_start+1);
        ofs << "" << "</Candidates>\n";
        //label_start++;

		//ofs << arg[j].putStringTypeOfSystematicAbsences() + " ("
			//	+ arg[j].putLatticeFigureOfMerit().putFiguresOfMerit().putLabel_FigureOfMeritWolff() + "="
				//+ num2str(arg[j].putLatticeFigureOfMerit().putFiguresOfMerit().putFigureOfMeritWolff(), 3) + ")\n";
	}
    /*ofs << "\n";
    label_start--;
    ofs.width(label_start);
    ofs << "" << "</TypeOfReflectionConditions>\n";*/

    ofs.width(label_start);
    ofs << "" << "</TypeOfReflectionConditions>\n";

    label_start--;
    ofs.width(label_start);
    ofs << "" << "</ConographOutput>\n";

    label_start--;
    ofs.width(label_start);
    ofs << "" << "</ZCodeParameters>\n";

	  ofs.close();
}


void printPeakPosition(
		const ControlParam& cdata,
		const PeakPosData& pdata, const Double& MIN_FOM,
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

#ifdef DEBUG
const DataReflectionConditions& ref_data = latfit.putDataOnReflectionConditions();
if( !(ref_data.isNotExtinct(hkl[0], hkl[1], hkl[2])) )
{
	cout << hkl[0] << " " << hkl[1] << " " << hkl[2] << "\n";
	cout << latfit.putStringTypeOfSystematicAbsences() << "\n";
	cout << latfit.putStringReflectionConditions() << "\n";
	assert(false);
}
#endif
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
		if( latfit_tray[j].putLatticeFigureOfMerit().putFiguresOfMerit().putFigureOfMeritWolff() < MIN_FOM ) break;
		*os << "X AppendToGraph yphase_" << j1 << " vs xphase_" << j1 << endl;
		*os << "X ModifyGraph offset(yphase_" << j1 << ")={0," << offset << "},mode(yphase_" << j1 << ")=3,marker(yphase_" << j1;
		*os << ")=10,msize(yphase_" << j1 << ")=3,mrkThick(yphase_" << j1 << ")=0.6,rgb(yphase_" << j1 << ")=(3,52428,1)" << endl;
	}

    *os << "X Legend/C/N=text0/J/A=MC \"";
    if( latfit_tray[0].putLatticeFigureOfMerit().putFiguresOfMerit().putFigureOfMeritWolff() >= MIN_FOM )
    {
    	*os << "\\s(yphase_" << 1 << ") "
    			+ latfit_tray[0].putShortStringTypeOfSystematicAbsences() + " ("
    			+ latfit_tray[0].putLatticeFigureOfMerit().putFiguresOfMerit().putLabel_FigureOfMeritWolff() + "="
    			+ num2str(latfit_tray[0].putLatticeFigureOfMerit().putFiguresOfMerit().putFigureOfMeritWolff(), 3) + ")";
    }
	for (Int4 j=1, j1=2; j<isize; j++, j1++)
	{
		if( latfit_tray[j].putLatticeFigureOfMerit().putFiguresOfMerit().putFigureOfMeritWolff() < MIN_FOM ) break;
		*os << "\\r\\s(yphase_" << j1 << ") "
				+ latfit_tray[j].putShortStringTypeOfSystematicAbsences() + " ("
				+ latfit_tray[j].putLatticeFigureOfMerit().putFiguresOfMerit().putLabel_FigureOfMeritWolff() + "="
				+ num2str(latfit_tray[j].putLatticeFigureOfMerit().putFiguresOfMerit().putFigureOfMeritWolff(), 3) + ")";
	}
	*os << "\"\nX Legend/C/N=text0/J/A=RT/X=0.00/Y=0.00\n";

    ofs.close();
}

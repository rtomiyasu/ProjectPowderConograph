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
#include <sstream>
#include "utility_data_structure/range.hh"
#include "utility_func/zstring.hh"
#include "utility_func/zmath.hh"
#include "zlog/zlog.hh"
#include "PeakPosData.hh"

using namespace std;

PeakPosData::PeakPosData() 
{
}

PeakPosData::~PeakPosData()
{
}

ZErrorMessageReadingFile PeakPosData::readFile(const string& filename)
{
	m_XWave.clear();
    m_YInt.clear();
    m_YErr.clear();
	PeakPosX.clear();
    PeakPosY.clear();
    PeakWidth.clear();
    toUseFlag.clear();

	ifstream ifs( filename.c_str() );

    if (!ifs){
   		return ZErrorMessageReadingFile(filename,
   				ZErrorMessage( ZErrorFileNotFound, __FILE__ , __LINE__, __FUNCTION__));
    }


    int count = 0;
    string s;

	try{
    	//  read headers while count <= 0.
		string header;
    	while( count <= 0 && getnewline(ifs, s) == ZErrorNoError )
    	{
        	if(is_blank(s)) continue;
	       	replace(s.begin(),s.end(),',',' ');
        	istringstream iss(s);
        	iss >> header;

        	if(count == 0)
        	{
        		if (header != "IGOR") throw ZErrorNotIGORFile;
        		count = -1;
        	}
        	else if (header.substr(0,5) == "WAVES")
        	{
        		iss >> m_XWave_title;
        		if( iss.fail() ) throw 1;
        		iss >> m_YInt_title;
        		if( iss.fail() ) throw 1;
        		iss >> m_YErr_title;
        		if( iss.fail() ) throw ZErrorErrorsAreNotContained;
        		
        		count = 1;
        	}
        	else
        	{
        		count = 0;
        		break;
        	}
    	}
		if( count <= 0 )
		{
			return ZErrorMessageReadingFile(filename,
					ZErrorMessage("There is no headers", __FILE__, __LINE__, __FUNCTION__));
		}


   		//read Profile Data.
    	Double t;
    	while (getnewline(ifs, s) == ZErrorNoError)
    	{
        	if(is_blank(s)) continue;
        	if(s.substr(0,5) == "BEGIN") continue;
       		if(s.substr(0,3) == "END" ) break;

        	istringstream iss(s);

    		iss >> t;
    		if( iss.fail() ) throw 1;
    		m_XWave.push_back(t);

    		iss >> t;
    		if( iss.fail() ) throw 1;
    		m_YInt.push_back(t);
    		
    		iss >> t;
			if( iss.fail() ) throw 1;
			m_YErr.push_back(t);
    	}
    
		
		bool flag = false;
    	while( getnewline(ifs, s) == ZErrorNoError )
    	{
        	if(is_blank(s)) continue;
	       	replace(s.begin(),s.end(),',',' ');
        	istringstream iss(s);
        	iss >> header;
        	
        	if (header.substr(0,5) == "WAVES")
        	{
        		flag = true;
        	}
    		break;
    	}
    	if( !flag )
    	{
        	return ZErrorMessageReadingFile(filename,
        			ZErrorMessage(ZErrorPeakPositionsAreNotContained, "There is no headers for peak positions", __FILE__, __LINE__, __FUNCTION__));
    	}

		while (getnewline(ifs, s) == ZErrorNoError)
		{
			if(is_blank(s)) continue;
			if(s.substr(0,5) == "BEGIN") continue;
			if(s.substr(0,3) == "END" ) break;

			istringstream iss(s);
			
			iss >> count;
			if( iss.fail() )
			{
				return ZErrorMessageReadingFile(filename,
						ZErrorMessage(ZErrorFileFormatBroken, __FILE__, __LINE__, __FUNCTION__));
			}

			iss >> t;
			if( iss.fail() )
			{
	        	return ZErrorMessageReadingFile(filename,
	        			ZErrorMessage(ZErrorFileFormatBroken, __FILE__, __LINE__, __FUNCTION__));
			}
			PeakPosX.push_back(t);

			iss >> t;
			if( iss.fail() )
			{
	        	return ZErrorMessageReadingFile(filename,
	        			ZErrorMessage(ZErrorFileFormatBroken, __FILE__, __LINE__, __FUNCTION__));
			}
			PeakPosY.push_back(t);

			iss >> t;
			if( iss.fail() )
			{
	        	return ZErrorMessageReadingFile(filename,
	        			ZErrorMessage(ZErrorFileFormatBroken, __FILE__, __LINE__, __FUNCTION__));
			}
			PeakWidth.push_back(t);

			iss >> flag;
			if( iss.fail() )
			{
	        	return ZErrorMessageReadingFile(filename,
	        			ZErrorMessage(ZErrorFileFormatBroken, __FILE__, __LINE__, __FUNCTION__));
			}
			toUseFlag.push_back(flag);

			if( *(PeakWidth.rbegin()) <= 0.0 )
			{
				*( PeakWidth.rbegin() ) = 0;
				*( toUseFlag.rbegin() ) = false;
ZLOG_WARN( "Peak position with non-positive FWHM is ignored : "
		+ num2str( *( PeakPosX.rbegin() ) ) + " "
		+ num2str( *( PeakPosY.rbegin() ) ) + " "
		+ num2str( *( PeakWidth.rbegin() ) ) + " "
		+ num2str( *( toUseFlag.rbegin() ) ) + "\n" );
			}
		}
    	if( PeakPosX.empty() )
    	{
        	return ZErrorMessageReadingFile(filename,
        			ZErrorMessage(ZErrorPeakPositionsAreNotContained, "There is no data for peak positions", __FILE__, __LINE__, __FUNCTION__));
    	}
	}
	catch(const ZErrorType& num)
	{
        if(num == ZErrorErrorsAreNotContained)
        {
        	return ZErrorMessageReadingFile(filename, ZErrorMessage(ZErrorErrorsAreNotContained, "The number of columns is less than "+ num2str(NumContents), __FILE__, __LINE__, __FUNCTION__));
        }
        else // num == ZErrorNotIGORFile
        {
        	return ZErrorMessageReadingFile(filename, ZErrorMessage(ZErrorNotIGORFile, "Not Igor File", __FILE__, __LINE__, __FUNCTION__));
        }
	}
	catch(const Int4& num)
	{
        if(num == 1)
        {
        	return ZErrorMessageReadingFile(filename,
        			ZErrorMessage("The number of columns is less than "+num2str<Int4>(NumContents), __FILE__, __LINE__, __FUNCTION__));
        }
        else // num == 3.
        {
        	return ZErrorMessageReadingFile(filename,
        			ZErrorMessage("Not Igor File", __FILE__, __LINE__, __FUNCTION__));
        }
	}

    return ZErrorMessageReadingFile();
}



void PeakPosData::printData(ostream *os) const
{
	*os << m_XWave_title << ", " << m_YInt_title << endl;
    *os << "BEGIN" << endl;
    os->setf(ios::right);
    os->setf(ios::uppercase);
    os->setf(ios::showpoint);

    os->precision(6);
    for (UInt4 i=0; i<m_XWave.size(); i++)
    {
        os->unsetf(ios::scientific);
        os->width(10);
        *os << m_XWave[i];
            
        os->setf(ios::scientific);
        os->width(15);
   	    *os << m_YInt[i];

//      os->width(15);
// 	    *os << m_YErr[i];
        *os << endl;
    }
    *os << "END" << endl;

    *os << "WAVES/O peak, peakpos, height, FWHM, Flag" << endl;
    *os << "BEGIN" << endl;
    for (UInt4 j = 0; j < PeakPosX.size(); j++)
    {
        // OUTPUT
        os->precision();
        os->width(5);
        *os << j + 1;
        
        os->width(15);
        *os << PeakPosX[j];
        os->width(15);
        *os << PeakPosY[j];
        os->width(15);
        *os << PeakWidth[j];

        os->precision();
        os->width(5);
        *os << toUseFlag[j] << endl;
    }
    *os << "END" << endl;
}


Double PeakPosData::putMaxPeakHeightOfFirst20() const
{
	Double ans = 0.0;
	const Int4 iend = min(20, (Int4)PeakPosX.size());
	if( iend <= 0 ) return ans;

	for(Int4 i=0; i<iend; i++)
	{
		ans = max(ans, PeakPosY[i]);
	}
	return ans;
}

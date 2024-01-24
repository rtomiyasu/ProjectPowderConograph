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
#ifndef _PeakPosData_h_
#define _PeakPosData_h_

// PeakPosData.hh
#include "RietveldAnalysisTypes.hh"
#include "zerror_type/error_out.hh"


using namespace std;

class PeakPosData
{
private:
    enum{ NumContents = 3 };	// x-phase, y-phase(, error).

	Vec_DP m_XWave;
	Vec_DP m_YInt;
	Vec_DP m_YErr;
	string m_XWave_title;
	string m_YInt_title;
	string m_YErr_title;

	Vec_DP PeakPosX;
    Vec_DP PeakPosY;
    Vec_DP PeakWidth;
    Vec_BOOL toUseFlag;

public:
    PeakPosData();
    ~PeakPosData();

    // For GUI developers.
    inline void setPeakPosXData(const Vec_DP& arg){ PeakPosX=arg; };
	inline void setPeakPosYData(const Vec_DP& arg){ PeakPosY=arg; };
	inline void setPeakWidthData(const Vec_DP& arg){ PeakWidth=arg; };
	inline void setToUseFlag(const Vec_BOOL& arg){ toUseFlag=arg; };

	inline void setXColumn(const Vec_DP& arg){ m_XWave=arg; };
	inline void setYIntColumn(const Vec_DP& arg){ m_YInt=arg; };
	inline void setYErrorColumn(const Vec_DP& arg){ m_YErr=arg; };

	inline void setXColumnTitle(const string& arg){ m_XWave_title = arg; };
	inline void setYIntColumnTitle(const string& arg){ m_YInt_title = arg; };
	inline void setYErrorColumnTitle(const string& arg){ m_YErr_title = arg; };

	inline const Vec_DP& putPeakPosXData() const { return PeakPosX; };
	inline const Vec_DP& putPeakPosYData() const { return PeakPosY; };
	inline const Vec_DP& putPeakWidthData() const { return PeakWidth; };
	inline const Vec_BOOL& putToUseFlag() const { return toUseFlag; };

	inline const Vec_DP& putXColumn() const { return m_XWave; };
	inline const Vec_DP& putYIntColumn() const { return m_YInt; };
	inline const Vec_DP& putYErrorColumn() const { return m_YErr;};

	inline const string& putXColumnTitle() const { return m_XWave_title; };
	inline const string& putYIntColumnTitle() const { return m_YInt_title; };
	inline const string& putYErrorColumnTitle() const { return m_YErr_title;};

	ZErrorMessageReadingFile readFile(const string& fname);
	void printData(ostream *os) const;
	Double putMaxPeakHeightOfFirst20() const;
};

#endif

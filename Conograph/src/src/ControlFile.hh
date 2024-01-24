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
#ifndef _ControlFile_h_
#define _ControlFile_h_

// ControlFile.hh
#include "utility_rw_param/I_ReadData.hh"

using namespace std;

class ControlFile : public I_ReadData
{
private:
    static const pair<RWParamProperty, RWParamData<string> > CntParamFname_Data;  // File name
    static const pair<RWParamProperty, RWParamData<string> > PeakDataFname_Data;  // File name
    static const pair<RWParamProperty, RWParamData<string> > OutFname_Data;  // File name

    string CntParamFname;  // File name
    string PeakDataFname;  // File name
    string OutFname;  // File name

public:
    ControlFile();
    virtual ~ControlFile();

    inline void setControlParamFileName(const string& arg){ CntParamFname = arg; };
    inline void settPeakDataFileName(const string& arg){ PeakDataFname = arg; };
    inline void setOutputFileName(const string& arg){ OutFname = arg; };

    inline const string& putControlParamFileName() const { return CntParamFname; };
    inline const string& putPeakDataFileName() const { return PeakDataFname; };
    inline const string& putOutputFileName() const { return OutFname; };
    
    const string& putTagLabel() const;
    void setData(const RWParamProperty& parent_prop,
			vector<RWParam_void>& tray);
};

#endif

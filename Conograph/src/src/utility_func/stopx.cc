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
#include "stopx.hh"
#include "../zerror_type/error_out.hh"

static int interrupted = 0;

void SetSignal(int p_signame)
{
  if (signal(p_signame, SigHandler) == SIG_ERR)
  {
    throw ZErrorMessage("The signal could not be set : "+num2str(p_signame), __FILE__, __LINE__, __FUNCTION__);
  }

  return;
}

void SigHandler(int p_signame)
{
	interrupted = 1;
}

int putInterruptionSignal()
{
	return interrupted;
}

// for GUI
#include <sstream>
#include <iomanip>
#include <cmath>
using namespace std;
namespace
{
	const char *functionName = "";
	int         progIndiv    = 0;
	int         progTotal    = 0;
}
void  interruptCansel() { interrupted = 2; }
void  interruptSkip  () { interrupted = 1; }
void  interruptClear () { interrupted = 0; }
string getProgress()
{
	ostringstream oss;
	oss <<setw(3)<< ::progTotal <<setw(4)<< ::progIndiv <<"%: "<< ::functionName;
	return oss.str();
}
void setProgress(const int indiv, const int total, const int span, const char *funcName)
{
	::progIndiv    = indiv==0 ? 0 : max(::progIndiv, indiv);
	::progTotal    = total + ::progIndiv * span/100;
	::functionName = funcName;
}
bool IS_CANSELED()
{
	return interrupted == 2;
}
void CHECK_INTERRUPTION() throw(const char*)
{
	if(IS_CANSELED()) { throw("abort"); }
}

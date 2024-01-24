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
#ifndef _ETYPE_ID_HH_
#define _ETYPE_ID_HH_

#include "../RietveldAnalysisTypes.hh"
#include <string>

typedef char PHSChar; // '0', '1', '2', 'U'.  

typedef enum {	//for switch
	_ZRietveldIDUnknown = -1,
	_ZRietveldIDFixed = 0,
	_ZRietveldIDVary = 1,
	_ZRietveldIDDepend = 2,
	_ZRietveldIDPhase = 3,
} etype_ID;

typedef struct{
	etype_ID ID;
	string phase_label;	// empty unless type = _ZRietveldIDPhase.
} ZRietveld_ID;

class ZErrorMessage;

ZErrorMessage checkFlag(const etype_ID& c, const bool flag = false);
PHSChar IDToPHSChar(const etype_ID& ID);
etype_ID PHSCharToID(const PHSChar& IDCharacter);

#endif

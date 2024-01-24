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
#include "etype_ID.hh"
#include "../zerror_type/error_mes.hh"

static const Char ZRietveldIDFixed = '0';
static const Char ZRietveldIDVary = '1';
static const Char ZRietveldIDDepend = '2';

ZErrorMessage checkFlag(const etype_ID& c, const bool flag)
{
	if( c==_ZRietveldIDFixed ) 	return ZErrorMessage();
	else if( c == _ZRietveldIDVary ) return ZErrorMessage();
	else if( flag && ( c == _ZRietveldIDDepend ) ) return ZErrorMessage();
	else{
        if(flag) return ZErrorMessage(ZErrorArgument, "INPUT 0-2", __FILE__, __LINE__, __FUNCTION__);
        else return ZErrorMessage(ZErrorArgument, "INPUT 0-1", __FILE__, __LINE__, __FUNCTION__);
	}
}


etype_ID PHSCharToID(const PHSChar& IDCharacter)
{
	switch (IDCharacter) {
		case ZRietveldIDFixed:{
			return _ZRietveldIDFixed;
		}break;
		case ZRietveldIDVary:{
			return _ZRietveldIDVary;
		}break;
		case ZRietveldIDDepend:{
			return _ZRietveldIDDepend;
		}break;
		default:{
		}break;
	}
	return _ZRietveldIDUnknown;
}


PHSChar IDToPHSChar(const etype_ID& ID)
{
	switch (ID) {
		case _ZRietveldIDFixed:{
			return ZRietveldIDFixed;
		}break;
		case _ZRietveldIDVary:{
			return ZRietveldIDVary;
		}break;
		case _ZRietveldIDDepend:{
			return ZRietveldIDDepend;
		}break;
		default:{
		}break;
	}
	return 'U';
}

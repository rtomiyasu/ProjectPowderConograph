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
#ifndef ZParawError_H_
#define ZParawError_H_

#include "../RietveldAnalysisTypes.hh"

class ZParawError 
{
public:
    double value; //Parameter's value
    double error; //Error value. (standard deviation)
	
	ZParawError(){ value = 0.0; error = 0.0; };
	ZParawError(const ZParawError& rhs){ value = rhs.value; error = rhs.error; };
	ZParawError(const Double& theValue){ value = theValue; error = 0.0; };
    ZParawError(const Double& theValue, const Double& theError){ value = theValue; error = theError; };
	
    virtual ~ZParawError(){};

    ZParawError& operator=(const ZParawError& rhs){ if(this != &rhs){ value = rhs.value; error = rhs.error; } return *this; };
    ZParawError& operator=(const Double& t){ value = t; error = 0.0; return *this; };
	inline ZParawError& operator*=(const Double& rhs){ value *= rhs; error *= rhs; return *this; };
};

inline ZParawError operator*(const ZParawError& lhs, Double& rhs)
{ 
	ZParawError ans = lhs;
	ans *= rhs;
	return ans; 
}

#endif /*ZParawError_H_*/

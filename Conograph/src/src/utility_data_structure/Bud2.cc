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
#include "Bud2.hh"

Bud2::Bud2()
{
	m_K[0] = -1;
	m_K[1] = -1;
	m_K[2] = -1;
	m_K[3] = -1;
}


Bud2::~Bud2()
{
}

void Bud2::setIndex(const Int4& K1, const Int4& K2, const Int4& K3, const Int4& K4)
{
#ifdef DEBUG
	Int4 count=0;
	if( K1 < 0 ) count++;
	if( K2 < 0 ) count++;
	if( K3 < 0 ) count++;
	if( K4 < 0 ) count++;
	assert( count < 2 );
#endif

	m_K[0] = K1;
	m_K[1] = K2;
	m_K[2] = K3;
	m_K[3] = K4;
}

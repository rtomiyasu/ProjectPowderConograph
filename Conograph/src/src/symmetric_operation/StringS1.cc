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
#include <cmath>
#include "StringS1.hh"
#include "../utility_func/zmath.hh"

StringS1::StringS1() : m_xyz_coef(0.0), m_s1(0.0)
{
}

StringS1::StringS1(const StringS1& rhs) : m_xyz_coef(rhs.m_xyz_coef), m_s1(rhs.m_s1)
{
}

StringS1::StringS1(const S1& s1) : 	m_xyz_coef(0.0), m_s1(s1)
{
}

StringS1::StringS1(const VecDat3<Double>& xyz_coef, const S1& s1) : m_xyz_coef(xyz_coef), m_s1(s1)
{
}

StringS1::~StringS1()
{
}

StringS1& StringS1::operator=(const StringS1& rhs)
{
	if(this != &rhs){
		m_xyz_coef = rhs.m_xyz_coef;
		m_s1 = rhs.m_s1;
	}
	return *this;
}

StringS1& StringS1::operator=(const S1& rhs)
{
	m_xyz_coef = 0.0;
	m_s1 = rhs;
	return *this;
}

const StringS1& StringS1::putX()
{
	static const VecDat3<Double> pos(1.0, 0.0, 0.0);
	static const StringS1 x(pos, 0.0);
	return x;
}

const StringS1& StringS1::putY()
{
	static const VecDat3<Double> pos(0.0, 1.0, 0.0);
	static const StringS1 y(pos, 0.0);
	return y;
}

const StringS1& StringS1::putZ()
{
	static const VecDat3<Double> pos(0.0, 0.0, 1.0);
	static const StringS1 z(pos, 0.0);
	return z;
}

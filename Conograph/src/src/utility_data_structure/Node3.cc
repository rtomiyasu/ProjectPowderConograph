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
#include <set>
#include "Node3.hh"
#include "../utility_data_structure/VCData.hh"
#include "../utility_data_structure/index_set.hh"
#include "../utility_func/transform_sym_matrix.hh"

Node3::Node3()
{
	m_detS = 0.0;
	m_detS_var = 0.0;
}

Node3::Node3(const Node3& rhs)
{
	for(Int4 i=0; i<NUM_K; i++) m_K[i] = rhs.m_K[i];
	for(Int4 i=0; i<NUM_SUM_K; i++) m_sum_K[i] = rhs.m_sum_K[i];
	m_detS = rhs.m_detS;
	m_detS_var = rhs.m_detS_var;
}

Node3::~Node3()
{
}

Node3& Node3::operator=(const Node3& rhs)
{
	if(this != &rhs)
	{
		for(Int4 i=0; i<NUM_K; i++) m_K[i] = rhs.m_K[i];
		for(Int4 i=0; i<NUM_SUM_K; i++) m_sum_K[i] = rhs.m_sum_K[i];
		m_detS = rhs.m_detS;
		m_detS_var = rhs.m_detS_var;
	}
	return *this;
}


bool Node3::setBud0(const Double& cv2,
		const VCData& K0,
		const VCData& K1,
		const VCData& K2,
		const VCData& K3,
		const VCData& K01,
		const VCData& K02,
		SymMat<VCData>& S)
{
	this->K(0) = K0;
	this->K(1) = K1;
	this->K(2) = K2;

	this->SumK(0, 1) = K01;
	this->SumK(0, 2) = K02;
	// |l1+l2|^2 = |l0|^2 + |l1|^2 + |l2|^2 + |l0+l1+l2|^2 - |l0+l1|^2 - |l0+l2|^2
	this->SumK(1, 2) = K0 + K1 + K2 + K3 - K01 - K02;
	
	S = calculateS3(*this);
	m_detS = Determinant3( chToDouble(S) );
	if( m_detS <= 0.0 ) return false;
	return true;
}

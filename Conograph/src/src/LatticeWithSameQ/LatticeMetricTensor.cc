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
#include "../utility_lattice_reduction/put_Selling_reduced_lattice.hh"
#include "../lattice_symmetry/gather_q_of_Ndim_lattice.hh"
#include "../lattice_symmetry/HKL_Q.hh"
#include "../laue_group/LaueGroup.hh"
#include "../zlog/zlog.hh"
#include "Choleskydcmp.hh"
#include "LatticeMetricTensor.hh"


LatticeMetricTensor::LatticeMetricTensor(const Int4& arg1, const Double& arg2)
	: m_S( SymMat<Double>(arg1), NRMat<Int4>(arg1+1,arg1) ),
	  m_determ_S(arg2)
{
	m_figures_of_merit = 0.0;
	m_num_lattice_found = 0;
}


LatticeMetricTensor::LatticeMetricTensor(const SymMatNplus1N_Double& S)
	: m_S( SymMat<Double>(S.first.size()), NRMat<Int4>(S.second.nrows(), S.second.ncols()) )
{
	assert( S.first.size() <= 4 );
	assert( S.first.size() == S.second.ncols() );
	assert( S.second.nrows() == S.second.ncols() + 1 );
	m_figures_of_merit = 0.0;
	m_num_lattice_found = 0;
	this->setLatticeConstants(S);
}



// This method assumes that S.first is Buerger-reduced and S.second * S.first * Transpose(S.second) is Selling-reduced.
void LatticeMetricTensor::setLatticeConstants(const SymMatNplus1N_Double& S)
{
	m_S = S;
	m_determ_S = Determinant(S.first);

	m_figures_of_merit = 0.0;
	m_num_lattice_found = 0;
}



void LatticeMetricTensor::setFigureOfMeritWu(const vector<Double>& qobs_tray,
		const vector< vector<HKL_Q>::const_iterator>& closest_qcal_tray)
{
	assert( qobs_tray.size() == closest_qcal_tray.size() );

	const size_t& num_qobs = qobs_tray.size();

	Double actually_disc = 0.0;
	for(size_t k=0; k<num_qobs; k++)
	{
		actually_disc += fabs( qobs_tray[k] - closest_qcal_tray[k]->Q() );
	}
	actually_disc /= num_qobs;

	const vector<HKL_Q>::const_iterator it_begin = closest_qcal_tray[0];
	const vector<HKL_Q>::const_iterator it_end = closest_qcal_tray[num_qobs-1] + 1;

	Double sum_diff = 0.0, diff;
	for(vector<HKL_Q>::const_iterator it=it_begin+1; it<it_end; it++)
	{
		diff = it->Q() - (it-1)->Q();
		sum_diff += diff * diff;
	}

	// Calculate the figure of merit by Wu.
	m_figures_of_merit = sum_diff / ( 4.0 * actually_disc * ( (it_end - 1)->Q() - it_begin->Q() ) );
}

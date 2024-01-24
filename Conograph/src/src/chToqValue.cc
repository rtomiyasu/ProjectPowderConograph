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
#include "chToqValue.hh"
#include "utility_func/zmath.hh"
#include "utility_data_structure/index_set.hh"


static ZErrorMessage chXDataToQData_TOF(const Vec_DP& posdata, const Vec_DP& pos_sterr_data,
		const Vec_BOOL& effective_flag,
		const Vec_DP& conv_param, vector<QData>& qdata, vector<Int4>& pos_qindex)
{
	Double d, Q, Q_err, dQ_dpos;
	try
	{
		if( conv_param.size() < 2 ) throw 4;

		const Int4 peak_num = posdata.size();

		Vec_DP dconv;
		for(UInt4 i=1; i<conv_param.size(); i++) dconv.push_back( conv_param[i]*i );

		qdata.clear();
		pos_qindex.clear();

		for(Int4 k=0; k<peak_num; k++)
		{
			if( !effective_flag[k] ) continue;
			if( pos_sterr_data[k] < 0.0 ) throw 5;

			if( !cal_polynomial_solution( conv_param, posdata[k], 
									(posdata[k] - conv_param[0]) / conv_param[1], d ) ) throw 2;
			if( d<=0.0 ) throw 3;

			Q = 1.0 / (d*d);

			dQ_dpos = 2.0 * Q / ( d * put_polynomial_value( dconv, d ) );
			Q_err = dQ_dpos * pos_sterr_data[k];

			qdata.push_back( QData(Q, Q_err * Q_err) );
			pos_qindex.push_back(k);
		}
	}
	catch(const int& num)
	{
		if( num == 2 ) return nerror_arg( "TOO MANY ITERATION FOR THE CALCULATION OF Q", __FILE__, __LINE__, __FUNCTION__ );
		else if( num == 3 ) return nerror_out_range( "THE SOLUTION FOR D IS OUT RANGE", d, __FILE__, __LINE__, __FUNCTION__ );
		else if( num == 4 ) return nerror_arg( "THE NUMBER OF CONVERSION PARMETERS IS LESS THAN 2", __FILE__, __LINE__, __FUNCTION__ );
		else if( num == 5 ) return nerror_arg( "A PEAK WIDTH IS A NEGATIVE NUMBER", __FILE__, __LINE__, __FUNCTION__ );
    }
	return ZErrorMessage();
}


ZErrorMessage chXDataToQData_AD(const Double& pos, const Double& pos_fwhm,
		const Double& peak_shift_param_rad, const Double& wlength, QData& qdata)
{
	static const Double INV_SQRT2LOG2 = 1.0 / sqrt( 8.0 * log(2.0) );
	static const Double RadDeg = PI() / 180.0; // = pi / 180.0.
	const Double pos_sterr = pos_fwhm * INV_SQRT2LOG2;

	Double w_d2, Q, Q_err, dQ_dtheta2;

	try
	{
		const Double theta2_rad = pos * RadDeg;
		const Double theta2_err_rad = pos_sterr * RadDeg;

			// w_d2 := wlength / 2d = sin(theta).
			w_d2 = sin( 0.5*(theta2_rad + peak_shift_param_rad) );
			if( w_d2 < 0.0 || 1.0 < w_d2 ) throw 3;

			Q = 2.0 * w_d2 / wlength;
			Q *= Q;

			// dw_d2/dtheta2 = 0.5*cos( 0.5*(theta2_rad + peak_shift_param_rad) )=0.5*sqrt(1-(w_d2)^2).
			// dQ / dw_d2 = 8*w_d2 / (wlength*wlength).
			dQ_dtheta2 = 4.0*w_d2*sqrt(1-w_d2*w_d2) / (wlength*wlength);
			Q_err = dQ_dtheta2 * theta2_err_rad;

			qdata.q = Q;
			qdata.q_var = Q_err * Q_err;
	}
	catch(const int& num)
	{
		return nerror_out_range("SIN(THETA) IS OUT RANGE", w_d2, __FILE__, __LINE__, __FUNCTION__ );
    }

	return ZErrorMessage();
}


// Solve pos = 2*theta + Z.
ZErrorMessage chXDataToQData_AD(const Vec_DP& posdata, const Vec_DP& pos_fwhm_data,
		const Vec_BOOL& effective_flag,
		const Vec_DP& peak_shift_param_rad, const Double& wlength,
		vector<QData>& pos_qdata, vector<Int4>& pos_qindex)
{
	try
	{
		if( peak_shift_param_rad.size() != 1 ) throw 4;
		if( wlength <= 0.0 ) throw 6;

		const Int4 peak_num = posdata.size();

		QData q2;

		pos_qdata.clear();
		pos_qindex.clear();

		ZErrorMessage zerr;

		for(Int4 k=0; k<peak_num; k++)
		{
			if( !effective_flag[k] ) continue;
			if( pos_fwhm_data[k] < 0.0 ) throw 5;

			zerr = chXDataToQData_AD(posdata[k], pos_fwhm_data[k],
					peak_shift_param_rad[0], wlength, q2);
			if( zerr.putErrorType() != ZErrorNoError ) return zerr;

			pos_qdata.push_back(q2);
			pos_qindex.push_back(k);
		}
	}
	catch(const int& num)
	{
		if( num == 4 ) return nerror_arg( "THE NUMBER OF PEAK SHIFT PARMETERS DOES NOT EQUAL 1", __FILE__, __LINE__, __FUNCTION__ );
		else if( num == 5 ) return nerror_arg( "A PEAK WIDTH IS A NEGATIVE NUMBER", __FILE__, __LINE__, __FUNCTION__ );
		else if( num == 6 ) return nerror_arg( "THE WAVE LENGTH IS A NEGATIVE NUMBER", __FILE__, __LINE__, __FUNCTION__ );
    }

	return ZErrorMessage();
}


// On output, Qdata is sorted into ascending order.
ZErrorMessage chToqValue(const ControlParam& cData,
		const PeakPosData& pData, vector<QData>& pos_qdata, vector<Int4>& pos_qindex)
{
	ZErrorMessage zerr;
	if( cData.IsAngleDispersion() )
	{
		zerr = chXDataToQData_AD(pData.putPeakPosXData(), pData.putPeakWidthData(), pData.putToUseFlag(),
				cData.putPeakShiftParamRadian(), cData.putWaveLength(), pos_qdata, pos_qindex);
		if( zerr.putErrorType() != ZErrorNoError ) return zerr;
	}
	else
	{
		zerr = chXDataToQData_TOF(pData.putPeakPosXData(), pData.putPeakWidthData(), pData.putToUseFlag(),
				cData.putConvParam(), pos_qdata, pos_qindex);
		if( zerr.putErrorType() != ZErrorNoError ) return zerr;
	}

	// Sort into ascending order.
	Vec_INT order;
	sort_order( pos_qdata, order );
	permute(order, pos_qindex.begin());

	return ZErrorMessage();
}

// Returns delta 2 theta = 2 arctan( (sin(pos2_rad) - 2*sin(pos_rad)) / (-cos(pos2_rad) + 2*cos(pos_rad)) ).
static Double fit_zero_shift_from_pair_of_peaks(const Double& sin_pos_rad,
		const Double& sin_pos2_rad, const Int4& k)
{
	return atan( (sin_pos2_rad - sin_pos_rad*k) / (-sqrt( 1.0 - sin_pos2_rad*sin_pos2_rad ) + sqrt( 1.0 - sin_pos_rad*sin_pos_rad )*k) ) * 2.0;
}


void fitZeroPointShift(
		const PeakPosData& pdata,
		const Int4& MAX_NUM_ITER,
		vector< pair<Int4, Int4> >& pair_of_q_index,
		vector<Double>& zero_point_shift_deg)
{
	static const Double RadDeg = PI() / 180.0; // = pi / 180.0.
	static const Double DegRad = 180.0 / PI(); // = 180.0 / pi.
	static const Int4 Max_Peak_Index = 30;

	const Vec_DP& posdata = pdata.putPeakPosXData();
	const Vec_BOOL& use_pos_flag = pdata.putToUseFlag();
	const Int4 num_iter = min( min(MAX_NUM_ITER, (Int4)posdata.size()), Max_Peak_Index);

	pair_of_q_index.clear();
	zero_point_shift_deg.clear();

	Vec_INT pos_order;
	Vec_DP sin_pos_tray;
	for(UInt4 i=0; i<posdata.size(); i++)
	{
		if( !use_pos_flag[i] ) continue;
		pos_order.push_back(i);
		sin_pos_tray.push_back(sin(0.5*posdata[i]*RadDeg));
	}
	if( sin_pos_tray.empty() ) return;
	if( (Int4)sin_pos_tray.size() > Max_Peak_Index ) sin_pos_tray.resize(Max_Peak_Index);

	Vec_INT order;
	sort_order(sin_pos_tray, order);
	permute(order, pos_order.begin());

	Double delta_2theta_rad;
	Int4 iend, ibegin, i;
	const Double Max_sin_pos = *(sin_pos_tray.rbegin());
	for(i=0; i<num_iter; i++)
	{
		const Int4 ntimes = ifloor(Max_sin_pos / sin_pos_tray[i]);

		for(Int4 k=2; k<=ntimes; k++)
		{
			const vector<Double>::iterator it2 = lower_bound( sin_pos_tray.begin(), sin_pos_tray.end(), sin_pos_tray[i]*k );

			const vector<Double>::iterator it_begin = max( sin_pos_tray.begin(), it2-1 );
			const vector<Double>::iterator it_end = min( sin_pos_tray.end(), it2+1 );

			ibegin = distance( sin_pos_tray.begin(), upper_bound( it_begin, it_end, sin_pos_tray[i]*1.95) );
			iend = distance( sin_pos_tray.begin(), lower_bound( it_begin, it_end, sin_pos_tray[i]*2.05) );

			for(Int4 q2_index=ibegin; q2_index<iend; q2_index++)
			{
				delta_2theta_rad = fit_zero_shift_from_pair_of_peaks( sin_pos_tray[i], sin_pos_tray[q2_index], k );

				pair_of_q_index.push_back( pair<Int4, Int4>( pos_order[i], pos_order[q2_index] ) );
				zero_point_shift_deg.push_back( delta_2theta_rad * DegRad );
			}
		}
	}
}

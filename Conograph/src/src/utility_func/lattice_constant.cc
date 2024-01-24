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
#include "lattice_constant.hh"
#include "transform_sym_matrix.hh"
#include "zmath.hh"


static const Double DegRad = 180.0 / PI();
static const Double RadDeg = PI() /180.0; // = pi / 180.0;

// lattice_constant a, b, c, alpha, beta, gamma(deg).
void calCoParameter(const VecDat3<Double>& length_axis, const VecDat3<Double>& angle_axis, SymMat<Double>& S)
{
	assert( S.size() == 3 );
	
	VecDat3<Double> cos_angle;
	for(Int4 k=0; k<3; k++){
		if( angle_axis[k] == 120.0 ) cos_angle[k] = -0.5;
		else if( angle_axis[k] == 90.0 ) cos_angle[k] = 0.0;
		else if( angle_axis[k] == 60.0 ) cos_angle[k] = 0.5;
		else cos_angle[k] = cos( RadDeg*angle_axis[k] );
	}

	const Double v = length_axis[0]*length_axis[1]*length_axis[2];
	const Double vol_sq = v * v * ( 1.0-cos_angle[0]*cos_angle[0]-cos_angle[1]*cos_angle[1]-cos_angle[2]*cos_angle[2]
																		+ 2.0*cos_angle[0]*cos_angle[1]*cos_angle[2] );
	const Double inv_vol_sq = 1.0/vol_sq;

	const VecDat3<Double> length_inn( 	length_axis[1]*length_axis[2],
										length_axis[0]*length_axis[2],
										length_axis[0]*length_axis[1]);

	// S(0,0)=A*, S(1,1)=B*, S(2,2)=C*. 
	for(Int4 k=0; k<3; k++){
		S(k,k)=length_inn[k]*length_inn[k]*inv_vol_sq*(1.0-cos_angle[k]*cos_angle[k]);
	}
	S(0,1)=length_inn[0]*length_inn[1]*inv_vol_sq*(cos_angle[0]*cos_angle[1]-cos_angle[2]);	// =F* 
	S(0,2)=length_inn[0]*length_inn[2]*inv_vol_sq*(cos_angle[0]*cos_angle[2]-cos_angle[1]);	// =E*
	S(1,2)=length_inn[1]*length_inn[2]*inv_vol_sq*(cos_angle[1]*cos_angle[2]-cos_angle[0]);	// =D*
}

// S_covar is the covariant matrix on S(A*,B*,C*,D*,E*,F*).
void changeCovariantMatrixStoLatticeConstant(const SymMat<Double>& S_covar, const VecDat3<Double>& length_axis,
const VecDat3<Double>& cos_angle, const VecDat3<Double>& sin_angle, SymMat<Double>& ans)
{
	static const Int4 pn_lat_const = 6;

	assert(S_covar.size() == pn_lat_const);
	assert(ans.size() == pn_lat_const);

	const Double aa=0.5*length_axis[0]*length_axis[0];
	const Double bb=0.5*length_axis[1]*length_axis[1];
	const Double cc=0.5*length_axis[2]*length_axis[2];
	const Double ab=length_axis[0]*length_axis[1];
	const Double ac=length_axis[0]*length_axis[2];
	const Double bc=length_axis[1]*length_axis[2];
	
	const Double cos2alpha=cos_angle[0]*cos_angle[0];
	const Double cos2beta=cos_angle[1]*cos_angle[1];
	const Double cos2gamma=cos_angle[2]*cos_angle[2];
	const Double cos_alpha_beta=cos_angle[0]*cos_angle[1];
	const Double cos_alpha_gamma=cos_angle[0]*cos_angle[2];
	const Double cos_beta_gamma=cos_angle[1]*cos_angle[2];
	
	// Jacobian matrix from S to lattice constants.
	NRMat<Double> Jac(pn_lat_const, pn_lat_const);
	
	Jac[0][0] = length_axis[0] * aa;					// = da / dS(0,0)
	Jac[0][1] = length_axis[0] * ab * cos_angle[2];		// = da / dS(0,1) 
	Jac[0][2] = length_axis[0] * ac * cos_angle[1];		// = da / dS(0,2)
	Jac[0][3] = length_axis[0] * bb * cos2gamma;		// = da / dS(1,1)
	Jac[0][4] = length_axis[0] * bc * cos_beta_gamma;	// = da / dS(1,2) 
	Jac[0][5] = length_axis[0] * cc * cos2beta;			// = da / dS(2,2)

	Jac[1][0] = length_axis[1] * aa * cos2gamma; 		// = db / dS(0,0)
	Jac[1][1] = length_axis[1] * ab * cos_angle[2];  	// = db / dS(0,1)
	Jac[1][2] = length_axis[1] * ac * cos_alpha_gamma; 	// = db / dS(0,2) 
	Jac[1][3] = length_axis[1] * bb;  					// = db / dS(1,1)
	Jac[1][4] = length_axis[1] * bc * cos_angle[0]; 	// = db / dS(1,2) 
	Jac[1][5] = length_axis[1] * cc * cos2alpha; 		// = db / dS(2,2)

	Jac[2][0] = length_axis[2] * aa * cos2beta; 		// = dc / dS(0,0)
	Jac[2][1] = length_axis[2] * ab * cos_alpha_beta;	// = dc / dS(0,1) 
	Jac[2][2] = length_axis[2] * ac * cos_angle[1]; 	// = dc / dS(0,2)
	Jac[2][3] = length_axis[2] * bb * cos2alpha; 		// = dc / dS(1,1)
	Jac[2][4] = length_axis[2] * bc * cos_angle[0]; 	// = dc / dS(1,2)
	Jac[2][5] = length_axis[2] * cc;					// = dc / dS(2,2)

	Jac[3][0] = -DegRad / sin_angle[0] * aa * ( 2.0 * cos_beta_gamma - cos_angle[0] * (cos2beta + cos2gamma) ); // = dalpha / dS(0,0)
	Jac[3][1] = -DegRad * sin_angle[0] * ab * cos_angle[1]; // = dalpha / dS(0,1) 
	Jac[3][2] = -DegRad * sin_angle[0] * ac * cos_angle[2]; // = dalpha / dS(0,2)
	Jac[3][3] = -DegRad * sin_angle[0] * bb * cos_angle[0]; // = dalpha / dS(1,1)
	Jac[3][4] = -DegRad * sin_angle[0] * bc;   				// = dalpha / dS(1,2)
	Jac[3][5] = -DegRad * sin_angle[0] * cc * cos_angle[0]; // = dalpha / dS(2,2)

	Jac[4][0] = -DegRad * sin_angle[1] * aa * cos_angle[1]; // = dbeta / dS(0,0)  
	Jac[4][1] = -DegRad * sin_angle[1] * ab * cos_angle[0]; // = dbeta / dS(0,1) 
	Jac[4][2] = -DegRad * sin_angle[1] * ac;  				// = dbeta / dS(0,2)
	Jac[4][3] = -DegRad / sin_angle[1] * bb * ( 2.0 * cos_alpha_gamma - cos_angle[1] * (cos2alpha + cos2gamma) ); // = dbeta / dS(1,1)
	Jac[4][4] = -DegRad * sin_angle[1] * bc * cos_angle[2]; // = dbeta / dS(1,2)  
	Jac[4][5] = -DegRad * sin_angle[1] * cc * cos_angle[1]; // = dbeta / dS(2,2)

	Jac[5][0] = -DegRad * sin_angle[2] * aa * cos_angle[2]; // = dgamma / dS(0,0) 
	Jac[5][1] = -DegRad * sin_angle[2] * ab; 				// = dgamma / dS(0,1) 
	Jac[5][2] = -DegRad * sin_angle[2] * ac * cos_angle[0]; // = dgamma / dS(0,2) 
	Jac[5][3] = -DegRad * sin_angle[2] * bb * cos_angle[2]; // = dgamma / dS(1,1) 
	Jac[5][4] = -DegRad * sin_angle[2] * bc * cos_angle[1]; // = dgamma / dS(1,2)   
	Jac[5][5] = -DegRad / sin_angle[2] * cc * ( 2.0 * cos_alpha_beta - cos_angle[2] * (cos2alpha + cos2beta) ); // = dgamma / dS(2,2) 

	// ans = Jac * S_covar * transpose(Jac).
	ans = transform_sym_matrix(Jac, S_covar);
}


static void calLatticeConstant(const SymMat<Double>& S, VecDat3<Double>& length_axis, VecDat3<Double>& angle_axis,
VecDat3<Double>& cos_angle, VecDat3<Double>& sin_angle)
{
	assert( S.size() >= 3 );

	// a*, b*, c*
	const VecDat3<Double> co_length( sqrt( S(0,0) ), sqrt( S(1,1) ), sqrt( S(2,2) ) );
	// alpha*, beta*, gamma*
	const VecDat3<Double> cos_co_angle( S(1,2)/(co_length[1]*co_length[2]), S(0,2)/(co_length[0]*co_length[2]), S(0,1)/(co_length[0]*co_length[1]) );
	const VecDat3<Double> sin_co_angle( 	sqrt(1.0 - cos_co_angle[0]*cos_co_angle[0]),
											sqrt(1.0 - cos_co_angle[1]*cos_co_angle[1]),
											sqrt(1.0 - cos_co_angle[2]*cos_co_angle[2])	);
	const Double inv_vol = co_length[0]*co_length[1]*co_length[2]
							* sqrt( 1.0-cos_co_angle[0]*cos_co_angle[0]-cos_co_angle[1]*cos_co_angle[1]-cos_co_angle[2]*cos_co_angle[2]
																+ 2.0*cos_co_angle[0]*cos_co_angle[1]*cos_co_angle[2] );
	const Double vol = 1.0/inv_vol;

	// alpha, beta, gamma
	cos_angle[0] = (cos_co_angle[1]*cos_co_angle[2]-cos_co_angle[0])/(sin_co_angle[1]*sin_co_angle[2]);
	cos_angle[1] = (cos_co_angle[0]*cos_co_angle[2]-cos_co_angle[1])/(sin_co_angle[0]*sin_co_angle[2]);
	cos_angle[2] = (cos_co_angle[0]*cos_co_angle[1]-cos_co_angle[2])/(sin_co_angle[0]*sin_co_angle[1]);
	sin_angle[0] = sqrt(1.0 - cos_angle[0]*cos_angle[0]);
	sin_angle[1] = sqrt(1.0 - cos_angle[1]*cos_angle[1]);
	sin_angle[2] = sqrt(1.0 - cos_angle[2]*cos_angle[2]);

	length_axis[0] = co_length[1]*co_length[2]*sin_co_angle[0]*vol;
	length_axis[1] = co_length[0]*co_length[2]*sin_co_angle[1]*vol;
	length_axis[2] = co_length[0]*co_length[1]*sin_co_angle[2]*vol;
	angle_axis[0] = atan2(sin_angle[0], cos_angle[0]) * DegRad;
	angle_axis[1] = atan2(sin_angle[1], cos_angle[1]) * DegRad;
	angle_axis[2] = atan2(sin_angle[2], cos_angle[2]) * DegRad;
}

void calLatticeConstant(const SymMat<Double>& S, VecDat3<Double>& length_axis, VecDat3<Double>& angle_axis)
{
	// alpha, beta, gamma
	VecDat3<Double> cos_angle, sin_angle;
	calLatticeConstant(S, length_axis, angle_axis, cos_angle, sin_angle);
}

void calLatticeConstant(const SymMatWCovar& S, VecDat3<Double>& length_axis, VecDat3<Double>& angle_axis,
SymMat<Double>& LatConst_covar)
{
	// alpha, beta, gamma
	VecDat3<Double> cos_angle, sin_angle;
	calLatticeConstant(S.ValMat(), length_axis, angle_axis, cos_angle, sin_angle);
	changeCovariantMatrixStoLatticeConstant(S.CovMat(), length_axis, cos_angle, sin_angle, LatConst_covar);
}

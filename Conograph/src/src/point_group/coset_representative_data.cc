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
#include<cmath>
#include<set>
#include"coset_representative_data.hh"
#include"point_gp_data.hh"
#include"../zerror_type/error_out.hh"

typedef struct{
	eGroupToMaxSubgp name;
	ePointGroup original_group;
	ePointGroup max_subgp;
	Int4 index;
	eSymmetricOperation e_gen_rep; // The generator of the representatives.
	bool normal_flag;	// If true, normal subgroup.
} group_to_max_subgp;

static const Int4 ISIZE_DATA = 430;
static const group_to_max_subgp group_to_max_subgp_data[ISIZE_DATA] = {
		{ Ci_to_C1,Ci,C1,2,Inv,1 },
		{ C2_X_to_C1,C2_X,C1,2,C2X,1 },
		{ C2_Y_to_C1,C2_Y,C1,2,C2Y,1 },
		{ C2_Z_to_C1,C2_Z,C1,2,C2Z,1 },
		{ C2D_X0_rho_to_C1,C2D_X0_rho,C1,2,C2primeX0_rho,1 },
		{ C2D_X1_rho_to_C1,C2D_X1_rho,C1,2,C2primeX1_rho,1 },
		{ C2D_Y0_rho_to_C1,C2D_Y0_rho,C1,2,C2primeY0_rho,1 },
		{ C2D_Y1_rho_to_C1,C2D_Y1_rho,C1,2,C2primeY1_rho,1 },
		{ C2D_Z0_to_C1,C2D_Z0,C1,2,C2primeZ0,1 },
		{ C2D_Z1_to_C1,C2D_Z1,C1,2,C2primeZ1,1 },
		{ C2D_X0_hex_to_C1,C2D_X0_hex,C1,2,C2primeX0_hex,1 },
		{ C2D_X1_hex_to_C1,C2D_X1_hex,C1,2,C2primeX1_hex,1 },
		{ C2D_Y0_hex_to_C1,C2D_Y0_hex,C1,2,C2primeY0_hex,1 },
		{ C2D_Y1_hex_to_C1,C2D_Y1_hex,C1,2,C2primeY1_hex,1 },
		{ Cs_X_to_C1,Cs_X,C1,2,SigmaX,1 },
		{ Cs_Y_to_C1,Cs_Y,C1,2,SigmaY,1 },
		{ Cs_Z_to_C1,Cs_Z,C1,2,SigmaZ,1 },
		{ CsD_X0_rho_to_C1,CsD_X0_rho,C1,2,SigmaDX0_rho,1 },
		{ CsD_X1_rho_to_C1,CsD_X1_rho,C1,2,SigmaDX1_rho,1 },
		{ CsD_Y0_rho_to_C1,CsD_Y0_rho,C1,2,SigmaDY0_rho,1 },
		{ CsD_Y1_rho_to_C1,CsD_Y1_rho,C1,2,SigmaDY1_rho,1 },
		{ CsD_Z0_to_C1,CsD_Z0,C1,2,SigmaDZ0,1 },
		{ CsD_Z1_to_C1,CsD_Z1,C1,2,SigmaDZ1,1 },
		{ CsD_X0_hex_to_C1,CsD_X0_hex,C1,2,SigmaDX0_hex,1 },
		{ CsD_X1_hex_to_C1,CsD_X1_hex,C1,2,SigmaDX1_hex,1 },
		{ CsD_Y0_hex_to_C1,CsD_Y0_hex,C1,2,SigmaDY0_hex,1 },
		{ CsD_Y1_hex_to_C1,CsD_Y1_hex,C1,2,SigmaDY1_hex,1 },
		{ C2h_X_to_Cs_X,C2h_X,Cs_X,2,Inv,1 },
		{ C2h_X_to_C2_X,C2h_X,C2_X,2,Inv,1 },
		{ C2h_X_to_Ci,C2h_X,Ci,2,SigmaX,1 },
		{ C2h_Y_to_Cs_Y,C2h_Y,Cs_Y,2,Inv,1 },
		{ C2h_Y_to_C2_Y,C2h_Y,C2_Y,2,Inv,1 },
		{ C2h_Y_to_Ci,C2h_Y,Ci,2,SigmaY,1 },
		{ C2h_Z_to_Cs_Z,C2h_Z,Cs_Z,2,Inv,1 },
		{ C2h_Z_to_C2_Z,C2h_Z,C2_Z,2,Inv,1 },
		{ C2h_Z_to_Ci,C2h_Z,Ci,2,SigmaZ,1 },
		{ C2hD_X0_rho_to_CsD_X0_rho,C2hD_X0_rho,CsD_X0_rho,2,Inv,1 },
		{ C2hD_X0_rho_to_C2D_X0_rho,C2hD_X0_rho,C2D_X0_rho,2,Inv,1 },
		{ C2hD_X0_rho_to_Ci,C2hD_X0_rho,Ci,2,SigmaDX0_rho,1 },
		{ C2hD_X1_rho_to_CsD_X1_rho,C2hD_X1_rho,CsD_X1_rho,2,Inv,1 },
		{ C2hD_X1_rho_to_C2D_X1_rho,C2hD_X1_rho,C2D_X1_rho,2,Inv,1 },
		{ C2hD_X1_rho_to_Ci,C2hD_X1_rho,Ci,2,SigmaDX1_rho,1 },
		{ C2hD_Y0_rho_to_CsD_Y0_rho,C2hD_Y0_rho,CsD_Y0_rho,2,Inv,1 },
		{ C2hD_Y0_rho_to_C2D_Y0_rho,C2hD_Y0_rho,C2D_Y0_rho,2,Inv,1 },
		{ C2hD_Y0_rho_to_Ci,C2hD_Y0_rho,Ci,2,SigmaDY0_rho,1 },
		{ C2hD_Y1_rho_to_CsD_Y1_rho,C2hD_Y1_rho,CsD_Y1_rho,2,Inv,1 },
		{ C2hD_Y1_rho_to_C2D_Y1_rho,C2hD_Y1_rho,C2D_Y1_rho,2,Inv,1 },
		{ C2hD_Y1_rho_to_Ci,C2hD_Y1_rho,Ci,2,SigmaDY1_rho,1 },
		{ C2hD_Z0_to_CsD_Z0,C2hD_Z0,CsD_Z0,2,Inv,1 },
		{ C2hD_Z0_to_C2D_Z0,C2hD_Z0,C2D_Z0,2,Inv,1 },
		{ C2hD_Z0_to_Ci,C2hD_Z0,Ci,2,SigmaDZ0,1 },
		{ C2hD_Z1_to_CsD_Z1,C2hD_Z1,CsD_Z1,2,Inv,1 },
		{ C2hD_Z1_to_C2D_Z1,C2hD_Z1,C2D_Z1,2,Inv,1 },
		{ C2hD_Z1_to_Ci,C2hD_Z1,Ci,2,SigmaDZ1,1 },
		{ C2hD_X0_hex_to_CsD_X0_hex,C2hD_X0_hex,CsD_X0_hex,2,Inv,1 },
		{ C2hD_X0_hex_to_C2D_X0_hex,C2hD_X0_hex,C2D_X0_hex,2,Inv,1 },
		{ C2hD_X0_hex_to_Ci,C2hD_X0_hex,Ci,2,SigmaDX0_hex,1 },
		{ C2hD_X1_hex_to_CsD_X1_hex,C2hD_X1_hex,CsD_X1_hex,2,Inv,1 },
		{ C2hD_X1_hex_to_C2D_X1_hex,C2hD_X1_hex,C2D_X1_hex,2,Inv,1 },
		{ C2hD_X1_hex_to_Ci,C2hD_X1_hex,Ci,2,SigmaDX1_hex,1 },
		{ C2hD_Y0_hex_to_CsD_Y0_hex,C2hD_Y0_hex,CsD_Y0_hex,2,Inv,1 },
		{ C2hD_Y0_hex_to_C2D_Y0_hex,C2hD_Y0_hex,C2D_Y0_hex,2,Inv,1 },
		{ C2hD_Y0_hex_to_Ci,C2hD_Y0_hex,Ci,2,SigmaDY0_hex,1 },
		{ C2hD_Y1_hex_to_CsD_Y1_hex,C2hD_Y1_hex,CsD_Y1_hex,2,Inv,1 },
		{ C2hD_Y1_hex_to_C2D_Y1_hex,C2hD_Y1_hex,C2D_Y1_hex,2,Inv,1 },
		{ C2hD_Y1_hex_to_Ci,C2hD_Y1_hex,Ci,2,SigmaDY1_hex,1 },
		{ D2_to_C2_Z,D2,C2_Z,2,C2X,1 },
		{ D2_to_C2_Y,D2,C2_Y,2,C2Z,1 },
		{ D2_to_C2_X,D2,C2_X,2,C2Y,1 },
		{ D2prime_X_rho_to_C2D_X1_rho,D2prime_X_rho,C2D_X1_rho,2,C2X,1 },
		{ D2prime_X_rho_to_C2D_X0_rho,D2prime_X_rho,C2D_X0_rho,2,C2X,1 },
		{ D2prime_X_rho_to_C2_X,D2prime_X_rho,C2_X,2,C2primeX0_rho,1 },
		{ D2prime_Y_rho_to_C2D_Y1_rho,D2prime_Y_rho,C2D_Y1_rho,2,C2Y,1 },
		{ D2prime_Y_rho_to_C2D_Y0_rho,D2prime_Y_rho,C2D_Y0_rho,2,C2Y,1 },
		{ D2prime_Y_rho_to_C2_Y,D2prime_Y_rho,C2_Y,2,C2primeY0_rho,1 },
		{ D2prime_Z_to_C2D_Z1,D2prime_Z,C2D_Z1,2,C2Z,1 },
		{ D2prime_Z_to_C2D_Z0,D2prime_Z,C2D_Z0,2,C2Z,1 },
		{ D2prime_Z_to_C2_Z,D2prime_Z,C2_Z,2,C2primeZ0,1 },
		{ D2prime_X_hex_to_C2D_X1_hex,D2prime_X_hex,C2D_X1_hex,2,C2Z,1 },
		{ D2prime_X_hex_to_C2D_X0_hex,D2prime_X_hex,C2D_X0_hex,2,C2Z,1 },
		{ D2prime_X_hex_to_C2_Z,D2prime_X_hex,C2_Z,2,C2primeX0_hex,1 },
		{ D2prime_Y_hex_to_C2D_Y1_hex,D2prime_Y_hex,C2D_Y1_hex,2,C2Z,1 },
		{ D2prime_Y_hex_to_C2D_Y0_hex,D2prime_Y_hex,C2D_Y0_hex,2,C2Z,1 },
		{ D2prime_Y_hex_to_C2_Z,D2prime_Y_hex,C2_Z,2,C2primeY0_hex,1 },
		{ C2v_X_to_Cs_Z,C2v_X,Cs_Z,2,C2X,1 },
		{ C2v_X_to_Cs_Y,C2v_X,Cs_Y,2,C2X,1 },
		{ C2v_X_to_C2_X,C2v_X,C2_X,2,SigmaY,1 },
		{ C2v_Y_to_Cs_Z,C2v_Y,Cs_Z,2,C2Y,1 },
		{ C2v_Y_to_Cs_X,C2v_Y,Cs_X,2,C2Y,1 },
		{ C2v_Y_to_C2_Y,C2v_Y,C2_Y,2,SigmaZ,1 },
		{ C2v_Z_to_Cs_Y,C2v_Z,Cs_Y,2,C2Z,1 },
		{ C2v_Z_to_Cs_X,C2v_Z,Cs_X,2,C2Z,1 },
		{ C2v_Z_to_C2_Z,C2v_Z,C2_Z,2,SigmaX,1 },
		{ C2vprime_X_rho_to_CsD_X1_rho,C2vprime_X_rho,CsD_X1_rho,2,C2X,1 },
		{ C2vprime_X_rho_to_CsD_X0_rho,C2vprime_X_rho,CsD_X0_rho,2,C2X,1 },
		{ C2vprime_X_rho_to_C2_X,C2vprime_X_rho,C2_X,2,SigmaDX0_rho,1 },
		{ C2vprime_Y_rho_to_CsD_Y1_rho,C2vprime_Y_rho,CsD_Y1_rho,2,C2Y,1 },
		{ C2vprime_Y_rho_to_CsD_Y0_rho,C2vprime_Y_rho,CsD_Y0_rho,2,C2Y,1 },
		{ C2vprime_Y_rho_to_C2_Y,C2vprime_Y_rho,C2_Y,2,SigmaDY0_rho,1 },
		{ C2vprime_Z_to_CsD_Z1,C2vprime_Z,CsD_Z1,2,C2Z,1 },
		{ C2vprime_Z_to_CsD_Z0,C2vprime_Z,CsD_Z0,2,C2Z,1 },
		{ C2vprime_Z_to_C2_Z,C2vprime_Z,C2_Z,2,SigmaDZ0,1 },
		{ C2vprime_X_hex_to_CsD_X1_hex,C2vprime_X_hex,CsD_X1_hex,2,C2Z,1 },
		{ C2vprime_X_hex_to_CsD_X0_hex,C2vprime_X_hex,CsD_X0_hex,2,C2Z,1 },
		{ C2vprime_X_hex_to_C2_Z,C2vprime_X_hex,C2_Z,2,SigmaDX0_hex,1 },
		{ C2vprime_Y_hex_to_CsD_Y1_hex,C2vprime_Y_hex,CsD_Y1_hex,2,C2Z,1 },
		{ C2vprime_Y_hex_to_CsD_Y0_hex,C2vprime_Y_hex,CsD_Y0_hex,2,C2Z,1 },
		{ C2vprime_Y_hex_to_C2_Z,C2vprime_Y_hex,C2_Z,2,SigmaDY0_hex,1 },
		{ C2vD_X0_rho_to_CsD_X1_rho,C2vD_X0_rho,CsD_X1_rho,2,SigmaX,1 },
		{ C2vD_X0_rho_to_Cs_X,C2vD_X0_rho,Cs_X,2,C2primeX0_rho,1 },
		{ C2vD_X0_rho_to_C2D_X0_rho,C2vD_X0_rho,C2D_X0_rho,2,SigmaX,1 },
		{ C2vD_X1_rho_to_CsD_X0_rho,C2vD_X1_rho,CsD_X0_rho,2,SigmaX,1 },
		{ C2vD_X1_rho_to_Cs_X,C2vD_X1_rho,Cs_X,2,C2primeX1_rho,1 },
		{ C2vD_X1_rho_to_C2D_X1_rho,C2vD_X1_rho,C2D_X1_rho,2,SigmaX,1 },
		{ C2vD_Y0_rho_to_CsD_Y1_rho,C2vD_Y0_rho,CsD_Y1_rho,2,SigmaY,1 },
		{ C2vD_Y0_rho_to_Cs_Y,C2vD_Y0_rho,Cs_Y,2,C2primeY0_rho,1 },
		{ C2vD_Y0_rho_to_C2D_Y0_rho,C2vD_Y0_rho,C2D_Y0_rho,2,SigmaY,1 },
		{ C2vD_Y1_rho_to_CsD_Y0_rho,C2vD_Y1_rho,CsD_Y0_rho,2,SigmaY,1 },
		{ C2vD_Y1_rho_to_Cs_Y,C2vD_Y1_rho,Cs_Y,2,C2primeY1_rho,1 },
		{ C2vD_Y1_rho_to_C2D_Y1_rho,C2vD_Y1_rho,C2D_Y1_rho,2,SigmaY,1 },
		{ C2vD_Z0_to_CsD_Z1,C2vD_Z0,CsD_Z1,2,SigmaZ,1 },
		{ C2vD_Z0_to_Cs_Z,C2vD_Z0,Cs_Z,2,C2primeZ0,1 },
		{ C2vD_Z0_to_C2D_Z0,C2vD_Z0,C2D_Z0,2,SigmaZ,1 },
		{ C2vD_Z1_to_CsD_Z0,C2vD_Z1,CsD_Z0,2,SigmaZ,1 },
		{ C2vD_Z1_to_Cs_Z,C2vD_Z1,Cs_Z,2,C2primeZ1,1 },
		{ C2vD_Z1_to_C2D_Z1,C2vD_Z1,C2D_Z1,2,SigmaZ,1 },
		{ C2vD_X0_hex_to_CsD_X1_hex,C2vD_X0_hex,CsD_X1_hex,2,SigmaZ,1 },
		{ C2vD_X0_hex_to_Cs_Z,C2vD_X0_hex,Cs_Z,2,C2primeX0_hex,1 },
		{ C2vD_X0_hex_to_C2D_X0_hex,C2vD_X0_hex,C2D_X0_hex,2,SigmaZ,1 },
		{ C2vD_X1_hex_to_CsD_X0_hex,C2vD_X1_hex,CsD_X0_hex,2,SigmaZ,1 },
		{ C2vD_X1_hex_to_Cs_Z,C2vD_X1_hex,Cs_Z,2,C2primeX1_hex,1 },
		{ C2vD_X1_hex_to_C2D_X1_hex,C2vD_X1_hex,C2D_X1_hex,2,SigmaZ,1 },
		{ C2vD_Y0_hex_to_CsD_Y1_hex,C2vD_Y0_hex,CsD_Y1_hex,2,SigmaZ,1 },
		{ C2vD_Y0_hex_to_Cs_Z,C2vD_Y0_hex,Cs_Z,2,C2primeY0_hex,1 },
		{ C2vD_Y0_hex_to_C2D_Y0_hex,C2vD_Y0_hex,C2D_Y0_hex,2,SigmaZ,1 },
		{ C2vD_Y1_hex_to_CsD_Y0_hex,C2vD_Y1_hex,CsD_Y0_hex,2,SigmaZ,1 },
		{ C2vD_Y1_hex_to_Cs_Z,C2vD_Y1_hex,Cs_Z,2,C2primeY1_hex,1 },
		{ C2vD_Y1_hex_to_C2D_Y1_hex,C2vD_Y1_hex,C2D_Y1_hex,2,SigmaZ,1 },
		{ D2h_to_C2v_Z,D2h,C2v_Z,2,Inv,1 },
		{ D2h_to_C2v_Y,D2h,C2v_Y,2,Inv,1 },
		{ D2h_to_C2v_X,D2h,C2v_X,2,Inv,1 },
		{ D2h_to_D2,D2h,D2,2,Inv,1 },
		{ D2h_to_C2h_Z,D2h,C2h_Z,2,SigmaX,1 },
		{ D2h_to_C2h_Y,D2h,C2h_Y,2,SigmaZ,1 },
		{ D2h_to_C2h_X,D2h,C2h_X,2,SigmaY,1 },
		{ D2hprime_X_rho_to_C2vD_X1_rho,D2hprime_X_rho,C2vD_X1_rho,2,Inv,1 },
		{ D2hprime_X_rho_to_C2vD_X0_rho,D2hprime_X_rho,C2vD_X0_rho,2,Inv,1 },
		{ D2hprime_X_rho_to_C2vprime_X_rho,D2hprime_X_rho,C2vprime_X_rho,2,Inv,1 },
		{ D2hprime_X_rho_to_D2prime_X_rho,D2hprime_X_rho,D2prime_X_rho,2,Inv,1 },
		{ D2hprime_X_rho_to_C2hD_X1_rho,D2hprime_X_rho,C2hD_X1_rho,2,SigmaX,1 },
		{ D2hprime_X_rho_to_C2hD_X0_rho,D2hprime_X_rho,C2hD_X0_rho,2,SigmaX,1 },
		{ D2hprime_X_rho_to_C2h_X,D2hprime_X_rho,C2h_X,2,SigmaDX0_rho,1 },
		{ D2hprime_Y_rho_to_C2vD_Y1_rho,D2hprime_Y_rho,C2vD_Y1_rho,2,Inv,1 },
		{ D2hprime_Y_rho_to_C2vD_Y0_rho,D2hprime_Y_rho,C2vD_Y0_rho,2,Inv,1 },
		{ D2hprime_Y_rho_to_C2vprime_Y_rho,D2hprime_Y_rho,C2vprime_Y_rho,2,Inv,1 },
		{ D2hprime_Y_rho_to_D2prime_Y_rho,D2hprime_Y_rho,D2prime_Y_rho,2,Inv,1 },
		{ D2hprime_Y_rho_to_C2hD_Y1_rho,D2hprime_Y_rho,C2hD_Y1_rho,2,SigmaY,1 },
		{ D2hprime_Y_rho_to_C2hD_Y0_rho,D2hprime_Y_rho,C2hD_Y0_rho,2,SigmaY,1 },
		{ D2hprime_Y_rho_to_C2h_Y,D2hprime_Y_rho,C2h_Y,2,SigmaDY0_rho,1 },
		{ D2hprime_Z_to_C2vD_Z1,D2hprime_Z,C2vD_Z1,2,Inv,1 },
		{ D2hprime_Z_to_C2vD_Z0,D2hprime_Z,C2vD_Z0,2,Inv,1 },
		{ D2hprime_Z_to_C2vprime_Z,D2hprime_Z,C2vprime_Z,2,Inv,1 },
		{ D2hprime_Z_to_D2prime_Z,D2hprime_Z,D2prime_Z,2,Inv,1 },
		{ D2hprime_Z_to_C2hD_Z1,D2hprime_Z,C2hD_Z1,2,SigmaZ,1 },
		{ D2hprime_Z_to_C2hD_Z0,D2hprime_Z,C2hD_Z0,2,SigmaZ,1 },
		{ D2hprime_Z_to_C2h_Z,D2hprime_Z,C2h_Z,2,SigmaDZ0,1 },
		{ D2hprime_X_hex_to_C2vD_X1_hex,D2hprime_X_hex,C2vD_X1_hex,2,Inv,1 },
		{ D2hprime_X_hex_to_C2vD_X0_hex,D2hprime_X_hex,C2vD_X0_hex,2,Inv,1 },
		{ D2hprime_X_hex_to_C2vprime_X_hex,D2hprime_X_hex,C2vprime_X_hex,2,Inv,1 },
		{ D2hprime_X_hex_to_D2prime_X_hex,D2hprime_X_hex,D2prime_X_hex,2,Inv,1 },
		{ D2hprime_X_hex_to_C2hD_X1_hex,D2hprime_X_hex,C2hD_X1_hex,2,SigmaZ,1 },
		{ D2hprime_X_hex_to_C2hD_X0_hex,D2hprime_X_hex,C2hD_X0_hex,2,SigmaZ,1 },
		{ D2hprime_X_hex_to_C2h_Z,D2hprime_X_hex,C2h_Z,2,SigmaDX0_hex,1 },
		{ D2hprime_Y_hex_to_C2vD_Y1_hex,D2hprime_Y_hex,C2vD_Y1_hex,2,Inv,1 },
		{ D2hprime_Y_hex_to_C2vD_Y0_hex,D2hprime_Y_hex,C2vD_Y0_hex,2,Inv,1 },
		{ D2hprime_Y_hex_to_C2vprime_Y_hex,D2hprime_Y_hex,C2vprime_Y_hex,2,Inv,1 },
		{ D2hprime_Y_hex_to_D2prime_Y_hex,D2hprime_Y_hex,D2prime_Y_hex,2,Inv,1 },
		{ D2hprime_Y_hex_to_C2hD_Y1_hex,D2hprime_Y_hex,C2hD_Y1_hex,2,SigmaZ,1 },
		{ D2hprime_Y_hex_to_C2hD_Y0_hex,D2hprime_Y_hex,C2hD_Y0_hex,2,SigmaZ,1 },
		{ D2hprime_Y_hex_to_C2h_Z,D2hprime_Y_hex,C2h_Z,2,SigmaDY0_hex,1 },
		{ C3_hex_to_C1,C3_hex,C1,3,C3plus_hex,1 },
		{ C31_rho_to_C1,C31_rho,C1,3,C31plus_rho,1 },
		{ C32_rho_to_C1,C32_rho,C1,3,C32plus_rho,1 },
		{ C33_rho_to_C1,C33_rho,C1,3,C33plus_rho,1 },
		{ C34_rho_to_C1,C34_rho,C1,3,C34plus_rho,1 },
		{ C3i_hex_to_C3_hex,C3i_hex,C3_hex,2,Inv,1 },
		{ C3i_hex_to_Ci,C3i_hex,Ci,3,C3plus_hex,1 },
		{ C31i_rho_to_C31_rho,C31i_rho,C31_rho,2,Inv,1 },
		{ C31i_rho_to_Ci,C31i_rho,Ci,3,C31plus_rho,1 },
		{ C32i_rho_to_C32_rho,C32i_rho,C32_rho,2,Inv,1 },
		{ C32i_rho_to_Ci,C32i_rho,Ci,3,C32plus_rho,1 },
		{ C33i_rho_to_C33_rho,C33i_rho,C33_rho,2,Inv,1 },
		{ C33i_rho_to_Ci,C33i_rho,Ci,3,C33plus_rho,1 },
		{ C34i_rho_to_C34_rho,C34i_rho,C34_rho,2,Inv,1 },
		{ C34i_rho_to_Ci,C34i_rho,Ci,3,C34plus_rho,1 },
		{ D3_0_hex_to_C3_hex,D3_0_hex,C3_hex,2,C2primeX0_hex,1 },
		{ D3_0_hex_to_C2D_Z0,D3_0_hex,C2D_Z0,3,C3plus_hex,0 },
		{ D3_0_hex_to_C2D_Y0_hex,D3_0_hex,C2D_Y0_hex,3,C3plus_hex,0 },
		{ D3_0_hex_to_C2D_X0_hex,D3_0_hex,C2D_X0_hex,3,C3plus_hex,0 },
		{ D3_1_hex_to_C3_hex,D3_1_hex,C3_hex,2,C2primeX1_hex,1 },
		{ D3_1_hex_to_C2D_Y1_hex,D3_1_hex,C2D_Y1_hex,3,C3plus_hex,0 },
		{ D3_1_hex_to_C2D_X1_hex,D3_1_hex,C2D_X1_hex,3,C3plus_hex,0 },
		{ D3_1_hex_to_C2D_Z1,D3_1_hex,C2D_Z1,3,C3plus_hex,0 },
		{ D31_rho_to_C31_rho,D31_rho,C31_rho,2,C2primeX0_rho,1 },
		{ D31_rho_to_C2D_Z0,D31_rho,C2D_Z0,3,C31plus_rho,0 },
		{ D31_rho_to_C2D_Y0_rho,D31_rho,C2D_Y0_rho,3,C31plus_rho,0 },
		{ D31_rho_to_C2D_X0_rho,D31_rho,C2D_X0_rho,3,C31plus_rho,0 },
		{ D32_rho_to_C32_rho,D32_rho,C32_rho,2,C2primeX1_rho,1 },
		{ D32_rho_to_C2D_Z0,D32_rho,C2D_Z0,3,C32plus_rho,0 },
		{ D32_rho_to_C2D_Y1_rho,D32_rho,C2D_Y1_rho,3,C32plus_rho,0 },
		{ D32_rho_to_C2D_X1_rho,D32_rho,C2D_X1_rho,3,C32plus_rho,0 },
		{ D33_rho_to_C33_rho,D33_rho,C33_rho,2,C2primeX0_rho,1 },
		{ D33_rho_to_C2D_Z1,D33_rho,C2D_Z1,3,C33plus_rho,0 },
		{ D33_rho_to_C2D_Y1_rho,D33_rho,C2D_Y1_rho,3,C33plus_rho,0 },
		{ D33_rho_to_C2D_X0_rho,D33_rho,C2D_X0_rho,3,C33plus_rho,0 },
		{ D34_rho_to_C34_rho,D34_rho,C34_rho,2,C2primeX1_rho,1 },
		{ D34_rho_to_C2D_Z1,D34_rho,C2D_Z1,3,C34plus_rho,0 },
		{ D34_rho_to_C2D_Y0_rho,D34_rho,C2D_Y0_rho,3,C34plus_rho,0 },
		{ D34_rho_to_C2D_X1_rho,D34_rho,C2D_X1_rho,3,C34plus_rho,0 },
		{ C3v_0_hex_to_C3_hex,C3v_0_hex,C3_hex,2,SigmaDX0_hex,1 },
		{ C3v_0_hex_to_CsD_Z0,C3v_0_hex,CsD_Z0,3,C3plus_hex,0 },
		{ C3v_0_hex_to_CsD_Y0_hex,C3v_0_hex,CsD_Y0_hex,3,C3plus_hex,0 },
		{ C3v_0_hex_to_CsD_X0_hex,C3v_0_hex,CsD_X0_hex,3,C3plus_hex,0 },
		{ C3v_1_hex_to_C3_hex,C3v_1_hex,C3_hex,2,SigmaDX1_hex,1 },
		{ C3v_1_hex_to_CsD_Y1_hex,C3v_1_hex,CsD_Y1_hex,3,C3plus_hex,0 },
		{ C3v_1_hex_to_CsD_X1_hex,C3v_1_hex,CsD_X1_hex,3,C3plus_hex,0 },
		{ C3v_1_hex_to_CsD_Z1,C3v_1_hex,CsD_Z1,3,C3plus_hex,0 },
		{ C31v_rho_to_C31_rho,C31v_rho,C31_rho,2,SigmaDX0_rho,1 },
		{ C31v_rho_to_CsD_Z0,C31v_rho,CsD_Z0,3,C31plus_rho,0 },
		{ C31v_rho_to_CsD_Y0_rho,C31v_rho,CsD_Y0_rho,3,C31plus_rho,0 },
		{ C31v_rho_to_CsD_X0_rho,C31v_rho,CsD_X0_rho,3,C31plus_rho,0 },
		{ C32v_rho_to_C32_rho,C32v_rho,C32_rho,2,SigmaDX1_rho,1 },
		{ C32v_rho_to_CsD_Z0,C32v_rho,CsD_Z0,3,C32plus_rho,0 },
		{ C32v_rho_to_CsD_Y1_rho,C32v_rho,CsD_Y1_rho,3,C32plus_rho,0 },
		{ C32v_rho_to_CsD_X1_rho,C32v_rho,CsD_X1_rho,3,C32plus_rho,0 },
		{ C33v_rho_to_C33_rho,C33v_rho,C33_rho,2,SigmaDX0_rho,1 },
		{ C33v_rho_to_CsD_Z1,C33v_rho,CsD_Z1,3,C33plus_rho,0 },
		{ C33v_rho_to_CsD_Y1_rho,C33v_rho,CsD_Y1_rho,3,C33plus_rho,0 },
		{ C33v_rho_to_CsD_X0_rho,C33v_rho,CsD_X0_rho,3,C33plus_rho,0 },
		{ C34v_rho_to_C34_rho,C34v_rho,C34_rho,2,SigmaDX1_rho,1 },
		{ C34v_rho_to_CsD_Z1,C34v_rho,CsD_Z1,3,C34plus_rho,0 },
		{ C34v_rho_to_CsD_Y0_rho,C34v_rho,CsD_Y0_rho,3,C34plus_rho,0 },
		{ C34v_rho_to_CsD_X1_rho,C34v_rho,CsD_X1_rho,3,C34plus_rho,0 },
		{ D3d_0_hex_to_C3v_0_hex,D3d_0_hex,C3v_0_hex,2,Inv,1 },
		{ D3d_0_hex_to_D3_0_hex,D3d_0_hex,D3_0_hex,2,Inv,1 },
		{ D3d_0_hex_to_C3i_hex,D3d_0_hex,C3i_hex,2,C2primeX0_hex,1 },
		{ D3d_0_hex_to_C2hD_Z0,D3d_0_hex,C2hD_Z0,3,C3plus_hex,0 },
		{ D3d_0_hex_to_C2hD_Y0_hex,D3d_0_hex,C2hD_Y0_hex,3,C3plus_hex,0 },
		{ D3d_0_hex_to_C2hD_X0_hex,D3d_0_hex,C2hD_X0_hex,3,C3plus_hex,0 },
		{ D3d_1_hex_to_C3v_1_hex,D3d_1_hex,C3v_1_hex,2,Inv,1 },
		{ D3d_1_hex_to_D3_1_hex,D3d_1_hex,D3_1_hex,2,Inv,1 },
		{ D3d_1_hex_to_C3i_hex,D3d_1_hex,C3i_hex,2,C2primeX1_hex,1 },
		{ D3d_1_hex_to_C2hD_Y1_hex,D3d_1_hex,C2hD_Y1_hex,3,C3plus_hex,0 },
		{ D3d_1_hex_to_C2hD_X1_hex,D3d_1_hex,C2hD_X1_hex,3,C3plus_hex,0 },
		{ D3d_1_hex_to_C2hD_Z1,D3d_1_hex,C2hD_Z1,3,C3plus_hex,0 },
		{ D31d_rho_to_C31v_rho,D31d_rho,C31v_rho,2,Inv,1 },
		{ D31d_rho_to_D31_rho,D31d_rho,D31_rho,2,Inv,1 },
		{ D31d_rho_to_C31i_rho,D31d_rho,C31i_rho,2,C2primeX0_rho,1 },
		{ D31d_rho_to_C2hD_Z0,D31d_rho,C2hD_Z0,3,C31plus_rho,0 },
		{ D31d_rho_to_C2hD_Y0_rho,D31d_rho,C2hD_Y0_rho,3,C31plus_rho,0 },
		{ D31d_rho_to_C2hD_X0_rho,D31d_rho,C2hD_X0_rho,3,C31plus_rho,0 },
		{ D32d_rho_to_C32v_rho,D32d_rho,C32v_rho,2,Inv,1 },
		{ D32d_rho_to_D32_rho,D32d_rho,D32_rho,2,Inv,1 },
		{ D32d_rho_to_C32i_rho,D32d_rho,C32i_rho,2,C2primeX1_rho,1 },
		{ D32d_rho_to_C2hD_Z0,D32d_rho,C2hD_Z0,3,C32plus_rho,0 },
		{ D32d_rho_to_C2hD_Y1_rho,D32d_rho,C2hD_Y1_rho,3,C32plus_rho,0 },
		{ D32d_rho_to_C2hD_X1_rho,D32d_rho,C2hD_X1_rho,3,C32plus_rho,0 },
		{ D33d_rho_to_C33v_rho,D33d_rho,C33v_rho,2,Inv,1 },
		{ D33d_rho_to_D33_rho,D33d_rho,D33_rho,2,Inv,1 },
		{ D33d_rho_to_C33i_rho,D33d_rho,C33i_rho,2,C2primeX0_rho,1 },
		{ D33d_rho_to_C2hD_Z1,D33d_rho,C2hD_Z1,3,C33plus_rho,0 },
		{ D33d_rho_to_C2hD_Y1_rho,D33d_rho,C2hD_Y1_rho,3,C33plus_rho,0 },
		{ D33d_rho_to_C2hD_X0_rho,D33d_rho,C2hD_X0_rho,3,C33plus_rho,0 },
		{ D34d_rho_to_C34v_rho,D34d_rho,C34v_rho,2,Inv,1 },
		{ D34d_rho_to_D34_rho,D34d_rho,D34_rho,2,Inv,1 },
		{ D34d_rho_to_C34i_rho,D34d_rho,C34i_rho,2,C2primeX1_rho,1 },
		{ D34d_rho_to_C2hD_Z1,D34d_rho,C2hD_Z1,3,C34plus_rho,0 },
		{ D34d_rho_to_C2hD_Y0_rho,D34d_rho,C2hD_Y0_rho,3,C34plus_rho,0 },
		{ D34d_rho_to_C2hD_X1_rho,D34d_rho,C2hD_X1_rho,3,C34plus_rho,0 },
		{ C4_X_to_C2_X,C4_X,C2_X,2,C4Xplus,1 },
		{ C4_Y_to_C2_Y,C4_Y,C2_Y,2,C4Yplus,1 },
		{ C4_Z_to_C2_Z,C4_Z,C2_Z,2,C4Zplus,1 },
		{ S4_X_to_C2_X,S4_X,C2_X,2,S4Xminus,1 },
		{ S4_Y_to_C2_Y,S4_Y,C2_Y,2,S4Yminus,1 },
		{ S4_Z_to_C2_Z,S4_Z,C2_Z,2,S4Zminus,1 },
		{ C4h_X_to_S4_X,C4h_X,S4_X,2,Inv,1 },
		{ C4h_X_to_C4_X,C4h_X,C4_X,2,Inv,1 },
		{ C4h_X_to_C2h_X,C4h_X,C2h_X,2,C4Xplus,1 },
		{ C4h_Y_to_S4_Y,C4h_Y,S4_Y,2,Inv,1 },
		{ C4h_Y_to_C4_Y,C4h_Y,C4_Y,2,Inv,1 },
		{ C4h_Y_to_C2h_Y,C4h_Y,C2h_Y,2,C4Yplus,1 },
		{ C4h_Z_to_S4_Z,C4h_Z,S4_Z,2,Inv,1 },
		{ C4h_Z_to_C4_Z,C4h_Z,C4_Z,2,Inv,1 },
		{ C4h_Z_to_C2h_Z,C4h_Z,C2h_Z,2,C4Zplus,1 },
		{ D4_X_to_C4_X,D4_X,C4_X,2,C2Y,1 },
		{ D4_X_to_D2prime_X_rho,D4_X,D2prime_X_rho,2,C4Xplus,1 },
		{ D4_X_to_D2,D4_X,D2,2,C4Xplus,1 },
		{ D4_Y_to_C4_Y,D4_Y,C4_Y,2,C2Z,1 },
		{ D4_Y_to_D2prime_Y_rho,D4_Y,D2prime_Y_rho,2,C4Yplus,1 },
		{ D4_Y_to_D2,D4_Y,D2,2,C4Yplus,1 },
		{ D4_Z_to_C4_Z,D4_Z,C4_Z,2,C2X,1 },
		{ D4_Z_to_D2prime_Z,D4_Z,D2prime_Z,2,C4Zplus,1 },
		{ D4_Z_to_D2,D4_Z,D2,2,C4Zplus,1 },
		{ C4v_X_to_C4_X,C4v_X,C4_X,2,SigmaY,1 },
		{ C4v_X_to_C2vprime_X_rho,C4v_X,C2vprime_X_rho,2,C4Xplus,1 },
		{ C4v_X_to_C2v_X,C4v_X,C2v_X,2,C4Xplus,1 },
		{ C4v_Y_to_C4_Y,C4v_Y,C4_Y,2,SigmaZ,1 },
		{ C4v_Y_to_C2vprime_Y_rho,C4v_Y,C2vprime_Y_rho,2,C4Yplus,1 },
		{ C4v_Y_to_C2v_Y,C4v_Y,C2v_Y,2,C4Yplus,1 },
		{ C4v_Z_to_C4_Z,C4v_Z,C4_Z,2,SigmaX,1 },
		{ C4v_Z_to_C2vprime_Z,C4v_Z,C2vprime_Z,2,C4Zplus,1 },
		{ C4v_Z_to_C2v_Z,C4v_Z,C2v_Z,2,C4Zplus,1 },
		{ D2d_X_to_S4_X,D2d_X,S4_X,2,C2Y,1 },
		{ D2d_X_to_C2vprime_X_rho,D2d_X,C2vprime_X_rho,2,S4Xminus,1 },
		{ D2d_X_to_D2,D2d_X,D2,2,S4Xminus,1 },
		{ D2d_Y_to_S4_Y,D2d_Y,S4_Y,2,C2Z,1 },
		{ D2d_Y_to_C2vprime_Y_rho,D2d_Y,C2vprime_Y_rho,2,S4Yminus,1 },
		{ D2d_Y_to_D2,D2d_Y,D2,2,S4Yminus,1 },
		{ D2d_Z_to_S4_Z,D2d_Z,S4_Z,2,C2X,1 },
		{ D2d_Z_to_C2vprime_Z,D2d_Z,C2vprime_Z,2,S4Zminus,1 },
		{ D2d_Z_to_D2,D2d_Z,D2,2,S4Zminus,1 },
		{ D2dprime_X_to_S4_X,D2dprime_X,S4_X,2,SigmaY,1 },
		{ D2dprime_X_to_C2v_X,D2dprime_X,C2v_X,2,S4Xminus,1 },
		{ D2dprime_X_to_D2prime_X_rho,D2dprime_X,D2prime_X_rho,2,S4Xminus,1 },
		{ D2dprime_Y_to_S4_Y,D2dprime_Y,S4_Y,2,SigmaZ,1 },
		{ D2dprime_Y_to_C2v_Y,D2dprime_Y,C2v_Y,2,S4Yminus,1 },
		{ D2dprime_Y_to_D2prime_Y_rho,D2dprime_Y,D2prime_Y_rho,2,S4Yminus,1 },
		{ D2dprime_Z_to_S4_Z,D2dprime_Z,S4_Z,2,SigmaX,1 },
		{ D2dprime_Z_to_C2v_Z,D2dprime_Z,C2v_Z,2,S4Zminus,1 },
		{ D2dprime_Z_to_D2prime_Z,D2dprime_Z,D2prime_Z,2,S4Zminus,1 },
		{ D4h_X_to_D2dprime_X,D4h_X,D2dprime_X,2,Inv,1 },
		{ D4h_X_to_D2d_X,D4h_X,D2d_X,2,Inv,1 },
		{ D4h_X_to_C4v_X,D4h_X,C4v_X,2,Inv,1 },
		{ D4h_X_to_D4_X,D4h_X,D4_X,2,Inv,1 },
		{ D4h_X_to_C4h_X,D4h_X,C4h_X,2,C2Y,1 },
		{ D4h_X_to_D2hprime_X_rho,D4h_X,D2hprime_X_rho,2,C4Xplus,1 },
		{ D4h_X_to_D2h,D4h_X,D2h,2,C4Xplus,1 },
		{ D4h_Y_to_D2dprime_Y,D4h_Y,D2dprime_Y,2,Inv,1 },
		{ D4h_Y_to_D2d_Y,D4h_Y,D2d_Y,2,Inv,1 },
		{ D4h_Y_to_C4v_Y,D4h_Y,C4v_Y,2,Inv,1 },
		{ D4h_Y_to_D4_Y,D4h_Y,D4_Y,2,Inv,1 },
		{ D4h_Y_to_C4h_Y,D4h_Y,C4h_Y,2,C2Z,1 },
		{ D4h_Y_to_D2hprime_Y_rho,D4h_Y,D2hprime_Y_rho,2,C4Yplus,1 },
		{ D4h_Y_to_D2h,D4h_Y,D2h,2,C4Yplus,1 },
		{ D4h_Z_to_D2dprime_Z,D4h_Z,D2dprime_Z,2,Inv,1 },
		{ D4h_Z_to_D2d_Z,D4h_Z,D2d_Z,2,Inv,1 },
		{ D4h_Z_to_C4v_Z,D4h_Z,C4v_Z,2,Inv,1 },
		{ D4h_Z_to_D4_Z,D4h_Z,D4_Z,2,Inv,1 },
		{ D4h_Z_to_C4h_Z,D4h_Z,C4h_Z,2,C2X,1 },
		{ D4h_Z_to_D2hprime_Z,D4h_Z,D2hprime_Z,2,C4Zplus,1 },
		{ D4h_Z_to_D2h,D4h_Z,D2h,2,C4Zplus,1 },
		{ C6_to_C3_hex,C6,C3_hex,2,C2Z,1 },
		{ C6_to_C2_Z,C6,C2_Z,3,C3plus_hex,1 },
		{ C3h_to_C3_hex,C3h,C3_hex,2,SigmaZ,1 },
		{ C3h_to_Cs_Z,C3h,Cs_Z,3,C3plus_hex,1 },
		{ C6h_to_C3h,C6h,C3h,2,Inv,1 },
		{ C6h_to_C6,C6h,C6,2,Inv,1 },
		{ C6h_to_C3i_hex,C6h,C3i_hex,2,C2Z,1 },
		{ C6h_to_C2h_Z,C6h,C2h_Z,3,C3plus_hex,1 },
		{ D6_to_C6,D6,C6,2,C2primeX0_hex,1 },
		{ D6_to_D3_1_hex,D6,D3_1_hex,2,C2Z,1 },
		{ D6_to_D3_0_hex,D6,D3_0_hex,2,C2Z,1 },
		{ D6_to_D2prime_Z,D6,D2prime_Z,3,C3plus_hex,0 },
		{ D6_to_D2prime_Y_hex,D6,D2prime_Y_hex,3,C3plus_hex,0 },
		{ D6_to_D2prime_X_hex,D6,D2prime_X_hex,3,C3plus_hex,0 },
		{ C6v_to_C6,C6v,C6,2,SigmaDX0_hex,1 },
		{ C6v_to_C3v_1_hex,C6v,C3v_1_hex,2,C2Z,1 },
		{ C6v_to_C3v_0_hex,C6v,C3v_0_hex,2,C2Z,1 },
		{ C6v_to_C2vprime_Z,C6v,C2vprime_Z,3,C3plus_hex,0 },
		{ C6v_to_C2vprime_Y_hex,C6v,C2vprime_Y_hex,3,C3plus_hex,0 },
		{ C6v_to_C2vprime_X_hex,C6v,C2vprime_X_hex,3,C3plus_hex,0 },
		{ D3h_0_hex_to_C3h,D3h_0_hex,C3h,2,C2primeX0_hex,1 },
		{ D3h_0_hex_to_C3v_1_hex,D3h_0_hex,C3v_1_hex,2,SigmaZ,1 },
		{ D3h_0_hex_to_D3_0_hex,D3h_0_hex,D3_0_hex,2,SigmaZ,1 },
		{ D3h_0_hex_to_C2vD_Z0,D3h_0_hex,C2vD_Z0,3,C3plus_hex,0 },
		{ D3h_0_hex_to_C2vD_Y0_hex,D3h_0_hex,C2vD_Y0_hex,3,C3plus_hex,0 },
		{ D3h_0_hex_to_C2vD_X0_hex,D3h_0_hex,C2vD_X0_hex,3,C3plus_hex,0 },
		{ D3h_1_hex_to_C3h,D3h_1_hex,C3h,2,C2primeX1_hex,1 },
		{ D3h_1_hex_to_C3v_0_hex,D3h_1_hex,C3v_0_hex,2,SigmaZ,1 },
		{ D3h_1_hex_to_D3_1_hex,D3h_1_hex,D3_1_hex,2,SigmaZ,1 },
		{ D3h_1_hex_to_C2vD_Z1,D3h_1_hex,C2vD_Z1,3,C3plus_hex,0 },
		{ D3h_1_hex_to_C2vD_Y1_hex,D3h_1_hex,C2vD_Y1_hex,3,C3plus_hex,0 },
		{ D3h_1_hex_to_C2vD_X1_hex,D3h_1_hex,C2vD_X1_hex,3,C3plus_hex,0 },
		{ D6h_to_D3h_1_hex,D6h,D3h_1_hex,2,Inv,1 },
		{ D6h_to_D3h_0_hex,D6h,D3h_0_hex,2,Inv,1 },
		{ D6h_to_C6v,D6h,C6v,2,Inv,1 },
		{ D6h_to_D6,D6h,D6,2,Inv,1 },
		{ D6h_to_C6h,D6h,C6h,2,C2primeX0_hex,1 },
		{ D6h_to_D3d_1_hex,D6h,D3d_1_hex,2,SigmaZ,1 },
		{ D6h_to_D3d_0_hex,D6h,D3d_0_hex,2,SigmaZ,1 },
		{ D6h_to_D2hprime_Z,D6h,D2hprime_Z,3,C3plus_hex,0 },
		{ D6h_to_D2hprime_Y_hex,D6h,D2hprime_Y_hex,3,C3plus_hex,0 },
		{ D6h_to_D2hprime_X_hex,D6h,D2hprime_X_hex,3,C3plus_hex,0 },
		{ T_to_D2,T,D2,3,C31plus_rho,1 },
		{ T_to_C34_rho,T,C34_rho,4,D2dummy,0 },
		{ T_to_C33_rho,T,C33_rho,4,D2dummy,0 },
		{ T_to_C32_rho,T,C32_rho,4,D2dummy,0 },
		{ T_to_C31_rho,T,C31_rho,4,D2dummy,0 },
		{ Th_to_T,Th,T,2,Inv,1 },
		{ Th_to_D2h,Th,D2h,3,C31plus_rho,1 },
		{ Th_to_C34i_rho,Th,C34i_rho,4,D2dummy,0 },
		{ Th_to_C33i_rho,Th,C33i_rho,4,D2dummy,0 },
		{ Th_to_C32i_rho,Th,C32i_rho,4,D2dummy,0 },
		{ Th_to_C31i_rho,Th,C31i_rho,4,D2dummy,0 },
		{ O_to_T,O,T,2,C4Xplus,1 },
		{ O_to_D4_Z,O,D4_Z,3,C31plus_rho,0 },
		{ O_to_D4_Y,O,D4_Y,3,C31plus_rho,0 },
		{ O_to_D4_X,O,D4_X,3,C31plus_rho,0 },
		{ O_to_D34_rho,O,D34_rho,4,D2dummy,0 },
		{ O_to_D33_rho,O,D33_rho,4,D2dummy,0 },
		{ O_to_D32_rho,O,D32_rho,4,D2dummy,0 },
		{ O_to_D31_rho,O,D31_rho,4,D2dummy,0 },
		{ Td_to_T,Td,T,2,S4Xminus,1 },
		{ Td_to_D2d_Z,Td,D2d_Z,3,C31plus_rho,0 },
		{ Td_to_D2d_Y,Td,D2d_Y,3,C31plus_rho,0 },
		{ Td_to_D2d_X,Td,D2d_X,3,C31plus_rho,0 },
		{ Td_to_C34v_rho,Td,C34v_rho,4,D2dummy,0 },
		{ Td_to_C33v_rho,Td,C33v_rho,4,D2dummy,0 },
		{ Td_to_C32v_rho,Td,C32v_rho,4,D2dummy,0 },
		{ Td_to_C31v_rho,Td,C31v_rho,4,D2dummy,0 },
		{ Oh_to_Td,Oh,Td,2,Inv,1 },
		{ Oh_to_O,Oh,O,2,Inv,1 },
		{ Oh_to_Th,Oh,Th,2,C4Xplus,1 },
		{ Oh_to_D4h_Z,Oh,D4h_Z,3,C31plus_rho,0 },
		{ Oh_to_D4h_Y,Oh,D4h_Y,3,C31plus_rho,0 },
		{ Oh_to_D4h_X,Oh,D4h_X,3,C31plus_rho,0 },
		{ Oh_to_D34d_rho,Oh,D34d_rho,4,D2dummy,0 },
		{ Oh_to_D33d_rho,Oh,D33d_rho,4,D2dummy,0 },
		{ Oh_to_D32d_rho,Oh,D32d_rho,4,D2dummy,0 },
		{ Oh_to_D31d_rho,Oh,D31d_rho,4,D2dummy,0 }
	};
	

inline const group_to_max_subgp& change_enum_to_data(const eGroupToMaxSubgp& num){

	if( num == dummy )
	{
		throw nerror_arg("dummy", __FILE__, __LINE__, __FUNCTION__);
	}
	return group_to_max_subgp_data[size_t(num)];
}

const ePointGroup& enumUpperGroup(const eGroupToMaxSubgp& num)
{
	return change_enum_to_data(num).original_group;
}

const ePointGroup& enumLowerGroup(const eGroupToMaxSubgp& num)
{
	return change_enum_to_data(num).max_subgp;
}

const Int4& Index(const eGroupToMaxSubgp& num)
{
	return change_enum_to_data(num).index;
}

// If ePointGroup = C1, returns dummy.
eGroupToMaxSubgp enumMaxNormalSubgroup(const ePointGroup& epg)
{
	if( epg != C1 ){
		for(Int4 k=0; k<ISIZE_DATA; k++)
			if( group_to_max_subgp_data[k].original_group == epg
				&& group_to_max_subgp_data[k].normal_flag ) return group_to_max_subgp_data[k].name;
	}
	return dummy;
}

void enumMaxSubgroup(const ePointGroup& epg, vector<eGroupToMaxSubgp>& tray)
{
	tray.clear();
	if( epg != C1 ){
		Int4 k;
		for(k=0; k<ISIZE_DATA; k++)
			if( group_to_max_subgp_data[k].original_group == epg ){
				tray.push_back( group_to_max_subgp_data[k].name );
				break;
		}
		k++;
		while( k<ISIZE_DATA && group_to_max_subgp_data[k].original_group == epg )
			tray.push_back( group_to_max_subgp_data[k++].name );
	}
}

void enumMinUpperGroup(const ePointGroup& epg, vector<eGroupToMaxSubgp>& tray)
{
	tray.clear();
	for(Int4 k=0; k<ISIZE_DATA; k++)
		if( group_to_max_subgp_data[k].max_subgp == epg )
				tray.push_back( group_to_max_subgp_data[k].name );
}

// Returns the representatives of the equivalent classes which do not include the unit element.
bool CosetRepresentativeMaxSubgp(const eGroupToMaxSubgp& num, Int4& index, eSymmetricOperation& symop)
{
	const group_to_max_subgp data = change_enum_to_data(num);
	index = data.index;
	symop = data.e_gen_rep;
	if( !data.normal_flag ) return false;
	return true;
}

// Returns true if and only if egp is a subgroup of egp2.
// ord1, ord2 are the order of the group egp, egp2 respectively.
static bool IsSubgroup(const ePointGroup& egp, const Int4& ord1, const ePointGroup& egp2, const Int4& ord2)
{
	if( ord2 % ord1 != 0 ) return false;
	if( egp == egp2 ) return true;
	if( ord1 == ord2 ) return false;
	if( egp == C1 ) return true;
	if( egp2 == C1 ) return false;

	vector<eGroupToMaxSubgp> tray;
	enumMaxSubgroup(egp2, tray);

	for(UInt4 j=0; j<tray.size(); j++)
		if( IsSubgroup( egp, ord1, enumLowerGroup(tray[j]), ord2/Index(tray[j]) ) ) return true;
	return false;
}


// Returns true if and only if egp is a subgroup of egp2.
bool IsSubgroup(const ePointGroup& egp, const ePointGroup& egp2)
{
	return IsSubgroup( egp, Order(egp), egp2, Order(egp2) );
}


// On output, epg is the group generated by subgp1 and subgp2.
bool generateGroup(const ePointGroup& subgp1, const ePointGroup& subgp2, ePointGroup& epg)
{
	if( IsSubgroup( subgp2, subgp1 ) ){
		epg = subgp1;
		return true;
	}
	if( IsSubgroup( subgp1, subgp2 ) ){
		epg = subgp2;
		return true;
	}

	vector<eGroupToMaxSubgp> tray, tray2;
	enumMinUpperGroup(subgp1, tray);
	ePointGroup epg2;

	for(UInt4 k=0; k<tray.size(); k++){
		if( !generateGroup( enumUpperGroup(tray[k]), subgp2, epg ) ) continue;
		bool flag = true;
		while(flag){
			enumMaxSubgroup(epg, tray2);
			flag = false;
			for(UInt4 j=0; j<tray2.size(); j++){
				epg2 = enumLowerGroup(tray2[j]);
				if( IsSubgroup( subgp1, epg2 ) && IsSubgroup( subgp2, epg2 ) ){
					epg = epg2;
					flag = true;
					break;
				}
			}
		}
		return true;
	}
	return false;
}


bool generateGroup(const vector<ePointGroup>& gp_vec, ePointGroup& egp)
{
	if( gp_vec.empty() ){
		egp = C1;
		return true;
	}
	
	set<ePointGroup> gp_set;
	const Int4 isize = gp_vec.size();
	for(Int4 j=0; j<isize; j++) gp_set.insert( gp_vec[j] );
	
	egp = *gp_set.rbegin();
	ePointGroup egp2;
	while( gp_set.size() > 1 ){
		gp_set.erase(*gp_set.rbegin());
		if( !generateGroup(egp, *gp_set.rbegin(), egp2) ) return false;
		gp_set.erase(*gp_set.rbegin());
		gp_set.insert(egp2);
		egp = *gp_set.rbegin();
	}
	return true;
}


bool generateGroup(const vector<SymmetricOperation>& symop, ePointGroup& egp)
{
	const Int4 isize = symop.size();
	vector<ePointGroup> gp_vec(isize);
	for(Int4 j=0; j<isize; j++){
		gp_vec[j] = generateGroup(symop[j]);
	}
	return generateGroup(gp_vec, egp);
}

void enumSubgroup(const ePointGroup& epg, set<ePointGroup>& esubpg)
{
	esubpg.clear();
	esubpg.insert(epg);
	if( epg == C1 ) return; 
	
	vector<eGroupToMaxSubgp> max_subgp;
	ePointGroup epg2;
	set<ePointGroup> esubpg2;
	enumMaxSubgroup(epg, max_subgp);
	for(vector<eGroupToMaxSubgp>::const_iterator it=max_subgp.begin(); it!=max_subgp.end(); it++)
	{
		epg2 = enumLowerGroup(*it);
		enumSubgroup(epg2, esubpg2);
		esubpg.insert(esubpg2.begin(), esubpg2.end());
	}
}


void enumSubgroup(const ePointGroup& epg, const ePointGroup& esubgp, set<ePointGroup>& esubpg)
{
	esubpg.clear();
	if( esubgp == epg || !IsSubgroup(esubgp, epg) ) return;
	esubpg.insert(epg);
	if( epg == C1 ) return; 
	
	vector<eGroupToMaxSubgp> max_subgp;
	ePointGroup epg2;
	set<ePointGroup> esubpg2;
	enumMaxSubgroup(epg, max_subgp);
	for(vector<eGroupToMaxSubgp>::const_iterator it=max_subgp.begin(); it!=max_subgp.end(); it++)
	{
		epg2 = enumLowerGroup(*it);
		enumSubgroup(epg2, esubgp, esubpg2);
		esubpg.insert(esubpg2.begin(), esubpg2.end());
	}
}



//void putGenerator(const ePointGroup& epg, vector<eSymmetricOperation>& soptray)
//{
//	soptray.clear();
//	if( epg == C1 ) return;
//	
//	Int4 index;
//	eSymmetricOperation sop;
//	ePointGroup epg2;
//	eGroupToMaxSubgp epg_subepg2;
//	
//	do{
//		epg_subepg2 = enumMaxNormalSubgroup(epg);
//		CosetRepresentativeMaxSubgp(epg_subepg2, index, sop);
//		soptray.push_back(sop);
//		epg2 = enumLowerGroup(epg_subepg2);
//	} while( epg2 != C1 );
//}

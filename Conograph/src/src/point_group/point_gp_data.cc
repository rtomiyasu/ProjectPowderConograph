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
#include"point_gp_data.hh"

// static const Int4 Max_Length = 7;

typedef struct{
	string schonflies_symbol;
	string HermannMauguin_symbol;
	Int4 order;
	ePointGroup epoint_gp;
	ePointGroup elaue_gp;
} point_gp;

static const point_gp& change_enum_to_data(const ePointGroup& epg)
{
	static const Int4 ISIZE = 136;
	static const point_gp point_gp_data[ISIZE] = {
		{ "C_1","1",1,C1,Ci },
		{ "C_i","-1",2,Ci,Ci },
		{ "C_2(a-axis)","2(a-axis)",2,C2_X,C2h_X },
		{ "C_2(b-axis)","2(b-axis)",2,C2_Y,C2h_Y },
		{ "C_2(c-axis)","2(c-axis)",2,C2_Z,C2h_Z },
		{ "C_2","2",2,C2D_X0_rho,C2hD_X0_rho },
		{ "C_2","2",2,C2D_X1_rho,C2hD_X1_rho },
		{ "C_2","2",2,C2D_Y0_rho,C2hD_Y0_rho },
		{ "C_2","2",2,C2D_Y1_rho,C2hD_Y1_rho },
		{ "C_2","2",2,C2D_Z0,C2hD_Z0 },
		{ "C_2","2",2,C2D_Z1,C2hD_Z1 },
		{ "C_2","2",2,C2D_X0_hex,C2hD_X0_hex },
		{ "C_2","2",2,C2D_X1_hex,C2hD_X1_hex },
		{ "C_2","2",2,C2D_Y0_hex,C2hD_Y0_hex },
		{ "C_2","2",2,C2D_Y1_hex,C2hD_Y1_hex },
		{ "C_s(a-axis)","m(a-axis)",2,Cs_X,C2h_X },
		{ "C_s(b-axis)","m(b-axis)",2,Cs_Y,C2h_Y },
		{ "C_s(c-axis)","m(c-axis)",2,Cs_Z,C2h_Z },
		{ "C_s","m",2,CsD_X0_rho,C2hD_X0_rho },
		{ "C_s","m",2,CsD_X1_rho,C2hD_X1_rho },
		{ "C_s","m",2,CsD_Y0_rho,C2hD_Y0_rho },
		{ "C_s","m",2,CsD_Y1_rho,C2hD_Y1_rho },
		{ "C_s","m",2,CsD_Z0,C2hD_Z0 },
		{ "C_s","m",2,CsD_Z1,C2hD_Z1 },
		{ "C_s","m",2,CsD_X0_hex,C2hD_X0_hex },
		{ "C_s","m",2,CsD_X1_hex,C2hD_X1_hex },
		{ "C_s","m",2,CsD_Y0_hex,C2hD_Y0_hex },
		{ "C_s","m",2,CsD_Y1_hex,C2hD_Y1_hex },
		{ "C_{2h}(a-axis)","2/m(a-axis)",4,C2h_X,C2h_X },
		{ "C_{2h}(b-axis)","2/m(b-axis)",4,C2h_Y,C2h_Y },
		{ "C_{2h}(c-axis)","2/m(c-axis)",4,C2h_Z,C2h_Z },
		{ "C_{2h}","2/m",4,C2hD_X0_rho,C2hD_X0_rho },
		{ "C_{2h}","2/m",4,C2hD_X1_rho,C2hD_X1_rho },
		{ "C_{2h}","2/m",4,C2hD_Y0_rho,C2hD_Y0_rho },
		{ "C_{2h}","2/m",4,C2hD_Y1_rho,C2hD_Y1_rho },
		{ "C_{2h}","2/m",4,C2hD_Z0,C2hD_Z0 },
		{ "C_{2h}","2/m",4,C2hD_Z1,C2hD_Z1 },
		{ "C_{2h}","2/m",4,C2hD_X0_hex,C2hD_X0_hex },
		{ "C_{2h}","2/m",4,C2hD_X1_hex,C2hD_X1_hex },
		{ "C_{2h}","2/m",4,C2hD_Y0_hex,C2hD_Y0_hex },
		{ "C_{2h}","2/m",4,C2hD_Y1_hex,C2hD_Y1_hex },
		{ "D_2","222",4,D2h,D2h },
		{ "D_2","222",4,D2hprime_X_rho,D2hprime_X_rho },
		{ "D_2","222",4,D2hprime_Y_rho },
		{ "D_2","222",4,D2hprime_Z,D2hprime_Z },
		{ "D_2","222",4,D2hprime_X_hex,D2hprime_X_hex },
		{ "D_2","222",4,D2hprime_Y_hex,D2hprime_Y_hex },
		{ "C_{2v}(a-axis)","2mm",4,C2v_X,D2h },
		{ "C_{2v}(b-axis)","m2m",4,C2v_Y,D2h },
		{ "C_{2v}(c-axis)","mm2",4,C2v_Z,D2h },
		{ "C_{2v}","mm2",4,C2vprime_X_rho,D2hprime_X_rho },
		{ "C_{2v}","mm2",4,C2vprime_Y_rho,D2hprime_Y_rho },
		{ "C_{2v}","mm2",4,C2vprime_Z,D2hprime_Z },
		{ "C_{2v}","mm2",4,C2vprime_X_hex,D2hprime_X_hex },
		{ "C_{2v}","mm2",4,C2vprime_Y_hex,D2hprime_Y_hex },
		{ "C_{2v}","mm2",4,C2vD_X0_rho,D2hprime_X_rho },
		{ "C_{2v}","mm2",4,C2vD_X1_rho,D2hprime_X_rho },
		{ "C_{2v}","mm2",4,C2vD_Y0_rho,D2hprime_Y_rho },
		{ "C_{2v}","mm2",4,C2vD_Y1_rho,D2hprime_Y_rho },
		{ "C_{2v}","mm2",4,C2vD_Z0,D2hprime_Z },
		{ "C_{2v}","mm2",4,C2vD_Z1,D2hprime_Z },
		{ "C_{2v}","mm2",4,C2vD_X0_hex,D2hprime_X_hex },
		{ "C_{2v}","mm2",4,C2vD_X1_hex,D2hprime_X_hex },
		{ "C_{2v}","mm2",4,C2vD_Y0_hex,D2hprime_Y_hex },
		{ "C_{2v}","mm2",4,C2vD_Y1_hex,D2hprime_Y_hex },
		{ "D_{2h}","mmm",8,D2h,D2h },
		{ "D_{2h}","mmm",8,D2hprime_X_rho,D2hprime_X_rho },
		{ "D_{2h}","mmm",8,D2hprime_Y_rho,D2hprime_Y_rho },
		{ "D_{2h}","mmm",8,D2hprime_Z,D2hprime_Z },
		{ "D_{2h}","mmm",8,D2hprime_X_hex },
		{ "D_{2h}","mmm",8,D2hprime_Y_hex },
		{ "C_3","3",3,C3_hex,C3i_hex },
		{ "C_3","3",3,C31_rho,C31i_rho },
		{ "C_3","3",3,C32_rho,C32i_rho },
		{ "C_3","3",3,C33_rho,C33i_rho },
		{ "C_3","3",3,C34_rho,C34i_rho },
		{ "C_{3i}","-3",6,C3i_hex,C3i_hex },
		{ "C_{3i}","-3",6,C31i_rho,C31i_rho },
		{ "C_{3i}","-3",6,C32i_rho,C32i_rho },
		{ "C_{3i}","-3",6,C33i_rho,C33i_rho },
		{ "C_{3i}","-3",6,C34i_rho,C34i_rho },
		{ "D_3","3.2",6,D3_0_hex,D3d_0_hex },
		{ "D_3","32.",6,D3_1_hex,D3d_1_hex },
		{ "D_3","32",6,D31_rho,D31d_rho },
		{ "D_3","32",6,D32_rho,D32d_rho },
		{ "D_3","32",6,D33_rho,D33d_rho },
		{ "D_3","32",6,D34_rho,D34d_rho },
		{ "C_{3v}","3.m",6,C3v_0_hex,D3d_0_hex },
		{ "C_{3v}","3m.",6,C3v_1_hex,D3d_1_hex },
		{ "C_{3v}","3m",6,C31v_rho,D31d_rho },
		{ "C_{3v}","3m",6,C32v_rho,D32d_rho },
		{ "C_{3v}","3m",6,C33v_rho,D33d_rho },
		{ "C_{3v}","3m",6,C34v_rho,D34d_rho },
		{ "D_{3d}","-3.m",12,D3d_0_hex,D3d_0_hex },
		{ "D_{3d}","-3m.",12,D3d_1_hex,D3d_1_hex },
		{ "D_{3d}","-3m",12,D31d_rho,D31d_rho },
		{ "D_{3d}","-3m",12,D32d_rho,D32d_rho },
		{ "D_{3d}","-3m",12,D33d_rho,D33d_rho },
		{ "D_{3d}","-3m",12,D34d_rho,D34d_rho },
		{ "C_4(a-axis)","4(a-axis)",4,C4_X,C4h_X },
		{ "C_4(b-axis)","4(b-axis)",4,C4_Y,C4h_Y },
		{ "C_4(c-axis)","4(c-axis)",4,C4_Z,C4h_Z },
		{ "S_4(a-axis)","-4(a-axis)",4,S4_X,C4h_X },
		{ "S_4(b-axis)","-4(b-axis)",4,S4_Y,C4h_Y },
		{ "S_4(c-axis)","-4(c-axis)",4,S4_Z,C4h_Z },
		{ "C_{4h}(a-axis)","4/m(a-axis)",8,C4h_X,C4h_X },
		{ "C_{4h}(b-axis)","4/m(b-axis)",8,C4h_Y,C4h_Y },
		{ "C_{4h}(c-axis)","4/m(c-axis)",8,C4h_Z,C4h_Z },
		{ "D_4(a-axis)","422",8,D4_X,D4h_X },
		{ "D_4(b-axis)","242",8,D4_Y,D4h_Y },
		{ "D_4(c-axis)","224",8,D4_Z,D4h_Z },
		{ "C_{4v}(a-axis)","4mm",8,C4v_X,D4h_X },
		{ "C_{4v}(b-axis)","m4m",8,C4v_Y,D4h_Y },
		{ "C_{4v}(c-axis)","mm4",8,C4v_Z,D4h_Z },
		{ "D_{2d}(a-axis)","-42m",8,D2d_X,D4h_X },
		{ "D_{2d}(b-axis)","-42m",8,D2d_Y,D4h_Y },
		{ "D_{2d}(c-axis)","-42m",8,D2d_Z,D4h_Z },
		{ "D_{2d}","-42m",8,D2dprime_X,D4h_X },
		{ "D_{2d}","-42m",8,D2dprime_Y,D4h_Y },
		{ "D_{2d}","-42m",8,D2dprime_Z,D4h_Z },
		{ "D_{4h}(a-axis)","4/mmm",16,D4h_X,D4h_X },
		{ "D_{4h}(b-axis)","m4/mm",16,D4h_Y,D4h_Y },
		{ "D_{4h}(c-axis)","mm4/m",16,D4h_Z,D4h_Z },
		{ "C_6","6",6,C6,C6h },
		{ "C_{3h}","-6",6,C3h,C6h },
		{ "C_{6h}","6/m",12,C6h,C6h },
		{ "D_{6}","622",12,D6,D6h },
		{ "C_{6v}","6mm",12,C6v,D6h },
		{ "D_{3h}","-6m2",12,D3h_0_hex,D6h },
		{ "D_{3h}","-62m",12,D3h_1_hex,D6h },
		{ "D_{6h}","6/mmm",24,D6h,D6h },
		{ "T","23",12,T,Th },
		{ "T_h","m-3",24,Th,Th },
		{ "O","432",24,O,Oh },
		{ "T_d","-43m",24,Td,Oh },
		{ "O_h","m-3m",48,Oh,Oh }
	};
	return point_gp_data[Int4(epg)];
}
	

const ePointGroup& enumLaueGroup(const ePointGroup& epg)
{
	return change_enum_to_data(epg).elaue_gp;
}

const string& Name(const ePointGroup& epg)
{
	return change_enum_to_data(epg).HermannMauguin_symbol;
}

const Int4& Order(const ePointGroup& epg)
{
	return change_enum_to_data(epg).order;
}

const string& SchonfliesSymbol(const ePointGroup& epg)
{
	return change_enum_to_data(epg).schonflies_symbol;
}


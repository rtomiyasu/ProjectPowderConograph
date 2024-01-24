#ifndef _gather_qcal2_HH_
#define _gather_qcal2_HH_
// set_additonal_Q.hh

#include "../lattice_symmetry/HKL_Q.hh"
#include "../utility_data_structure/SymMat.hh"
//#include "../lattice_symmetry/gather_q_of_Ndim_lattice.hh"

using namespace std;

void gatherQcal(const SymMat<Double>& S_super,
		        const Double& maxQ,
				const NRMat<Int4>& transform_hkl,
                const BravaisType& ebrav_type,
                const Int4& irc_type,
		        vector<HKL_Q>& qcal_tray);

/*(const SymMat<Double>& S_super,
		const Double& maxQ, const eTypeOfSystematicAbsence& etype,
		vector<HKL_Q>& qcal_tray);*/


#endif

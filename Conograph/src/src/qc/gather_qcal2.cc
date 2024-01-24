#include "reflection_conditions.hh"
#include "gather_qcal2.hh"
#include "../lattice_symmetry/gather_q_of_Ndim_lattice.hh"

void gatherQcal(const SymMat<Double>& S_super,
		        const Double& maxQ,
				const NRMat<Int4>& transform_hkl,
                const BravaisType& brav_type,
                const Int4& irc_type,
		        vector<HKL_Q>& qcal_tray)
{
    // First call gatherQcal
	gatherQcal(S_super, maxQ, transform_hkl, qcal_tray);
	if( irc_type < 0 ) return;

    // Next, erase entries of qcal_tray, according to ebrav_type and erc_type.
    const DataReflectionConditions& data = putInformationOnReflectionConditions(brav_type, irc_type);
    Int4 index = 0;
    for(vector<HKL_Q>::const_iterator it=qcal_tray.begin(); it!=qcal_tray.end(); it++)
    {
    	const VecDat3<Int4>& hkl = it->HKL();
    	if( (data.isNotExtinct(hkl[0], hkl[1], hkl[2]) ) )
    	{
            	qcal_tray[index++] = *it;
    	}
    }
    qcal_tray.erase(qcal_tray.begin()+index, qcal_tray.end());
}

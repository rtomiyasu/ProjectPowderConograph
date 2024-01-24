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
#ifndef PGNORMALSERIESTRAY_HH_
#define PGNORMALSERIESTRAY_HH_

#include"../RietveldAnalysisTypes.hh"
#include"../point_group/enumPointGroup.hh"
#include"../symmetric_operation/MillerIndex.hh"
#include"../symmetric_operation/SymmetricOperation.hh"

template <class T>
class XYZCoord;
class S1;
class StringS1;

// Class of a normal series of a space-group.
class PGNormalSeriesTray
{
private:
	ePointGroup m_epg;
	
	Int4 HKLEquiv(const MillerIndex3&, const Int4&, 
						const vector<Int4>::const_iterator&,
						const vector<SymmetricOperation>::const_iterator&,
						vector<MillerIndex3>&) const;

protected:
	// Normal series of m_esg.
	Int4 m_length;      // The length of the normal series m_esg.
	
	vector<Int4> m_index; 
	vector< eSymmetricOperation > m_coset_rep;  // Coset representatives of the normal series(The size equal m_length - 1).
	
	inline vector<SymmetricOperation> putSymmetricOperation() const;

public:
	PGNormalSeriesTray(const ePointGroup&);
	~PGNormalSeriesTray();

	inline const ePointGroup& enumPointGroup() const;
	inline Int4 putHKLEquiv(const MillerIndex3&, vector<MillerIndex3>&) const;

};

inline const ePointGroup& PGNormalSeriesTray::enumPointGroup() const
{
	return m_epg;
}

inline vector<SymmetricOperation> PGNormalSeriesTray::putSymmetricOperation() const
{
	vector<SymmetricOperation> ans(m_length-1);
	for(Int4 k=0; k<m_length-1; k++)
	{
		ans[k] = change_enum_to_data(m_coset_rep[k]);
	}
	return ans;
}

inline Int4 PGNormalSeriesTray::putHKLEquiv(const MillerIndex3& hkl,
vector<MillerIndex3>& hkl_equiv) const
{
	return HKLEquiv(hkl, m_length, m_index.begin(), putSymmetricOperation().begin(), hkl_equiv);
}


#endif /*PGNORMALSERIESTRAY_HH_*/

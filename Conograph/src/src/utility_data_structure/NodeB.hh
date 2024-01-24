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
#ifndef _NodeB_hh_
#define _NodeB_hh_
// Node3.hh

#include <set>
#include "Bud.hh"
#include "VCData.hh"
#include "VecDat3.hh"
#include "../RietveldAnalysisTypes.hh"
#include "../utility_func/zstring.hh"


class NodeB
{
private:
	Int4 m_left;	// K1
	Int4 m_right;	// K2
	NodeB* m_left_branch;	// K1, K4
	NodeB* m_right_branch;	// K2, K4 
	
public:
	NodeB();
	NodeB(const Int4& i, const Int4& j);
	NodeB(const NodeB& i, const NodeB& j);
	NodeB(const NodeB& rhs);
    ~NodeB();

    NodeB& operator=(const NodeB& rhs);
    void set(NodeB& rhs);

    inline const Int4& Left() const;
    inline const Int4& Right() const;
    inline bool IsBud() const;
    inline Int4 Upper() const;
    
    void count(std::set<Int4>& index) const;

    inline const NodeB& LeftBranch() const;
    inline const NodeB& RightBranch() const;

    void setBranch(const Int4&);
    void cutBranch();

    void putRootBud(const Int4& K3, std::set<Bud>& tray) const;
    void putBud(const Int4& K3, std::set<Bud>& tray) const;
//    void putQuadraticForm(const VecDat3<Int4>& hkl_left,
//    						const VecDat3<Int4>& hkl_right,
//    						multimap<Int4, VecDat3<Int4> >& qindex_hkl) const;

    inline Int4 count_bud() const;

    inline string chToString() const;
    inline void swapBranch();
    void print(ostream&, const Int4& num, const Int4& K3, const Double& maxQ) const;
};


inline Int4 NodeB::count_bud() const
{
	if( IsBud() ) return 0;
	
	if( m_left == m_right ) return m_left_branch->count_bud() + 1;
	else return m_left_branch->count_bud() + m_right_branch->count_bud() + 1;
}


inline const Int4& NodeB::Left() const
{
	return m_left;
}


inline const Int4& NodeB::Right() const
{
	return m_right;
}


inline bool NodeB::IsBud() const
{
	return m_left_branch == NULL && m_right_branch == NULL;
}


inline Int4 NodeB::Upper() const
{
	if( this->IsBud() )
	{
   		throw ZErrorMessage(ZErrorNullPointer, __FILE__, __LINE__, __FUNCTION__);
	}
	return m_left_branch->Right();
}


inline const NodeB& NodeB::LeftBranch() const
{
	if(m_left_branch==NULL)
	{
   		throw ZErrorMessage(ZErrorNullPointer, __FILE__, __LINE__, __FUNCTION__);
	}
	return *m_left_branch;
}


inline const NodeB& NodeB::RightBranch() const
{
	if(m_right_branch==NULL)
	{
       	throw ZErrorMessage(ZErrorNullPointer, __FILE__, __LINE__, __FUNCTION__);
	}
	return *m_right_branch;
}


inline string NodeB::chToString() const
{
	if( m_left < 0 ) return "(-1, " + num2str<Int4>(m_right+1) + ")";
	if( m_right < 0 ) return "(" + num2str<Int4>(m_left+1) + ", -1)";
	return "(" + num2str<Int4>(m_left+1) + ", " + num2str<Int4>(m_right+1) + ")";
}

inline void NodeB::swapBranch()
{
	if( m_left == m_right ) return;
	swap(m_left, m_right);
	swap(m_left_branch, m_right_branch);
}


#endif

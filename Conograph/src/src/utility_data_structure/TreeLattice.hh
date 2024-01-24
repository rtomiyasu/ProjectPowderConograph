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
#ifndef _TreeLattice_hh_
#define _TreeLattice_hh_
// TreeLattice.hh

#include <ostream>
#include <fstream>
#include <assert.h>

#include "../RietveldAnalysisTypes.hh"
#include "NodeB.hh"
#include "VecDat3.hh"

class TreeLattice
{
private:
//	const Int4 m_data_type;

	NodeB* m_root;	// K1, K2
	NodeB* m_root_on_left;	// K2, K3
	NodeB* m_root_on_right;	// K3, K1
	
	bool HeadIsTail;

	bool m_is_set_sort_criteria;
	Int4 m_count_Q;
	Double m_detS;
	
    // Returns 2 if HeadIsTail is true.
    // Returns 1 if HeadIsTail is false, and m_root, m_root_on_left, m_root_on_right are not NULL.
    // Returns 0 otherwise.
//	inline Int4 SuperBasis() const;

	// Sets m_detS.
    void setAreaSquare();
	// Sets m_count_Q.
    void setCountOfQ();

public:
	TreeLattice();
	TreeLattice(const TreeLattice& rhs);
    ~TreeLattice();

    inline const NodeB& Root() const;

    inline void setRoot(const Int4&, const Int4&);
    inline void setRoot(const NodeB&, const NodeB&);

    inline void setRootEqualUpper();
    inline void setRootOnLeftBranch(const NodeB& node);
    inline void setRootOnRightBranch(const NodeB& node);

    inline void setSortingCriteria() { this->setCountOfQ(); this->setAreaSquare(); m_is_set_sort_criteria = true; };
    inline Int4 putCountOfQ() const{ assert(m_is_set_sort_criteria); return m_count_Q; };
    inline const Double& putAreaSquare() const{ assert(m_is_set_sort_criteria); return m_detS; };

    inline void swapBranch();

    TreeLattice& operator=(const TreeLattice& rhs);
   
    void clear();

    // On output, index_tray is sorted into ascending order, any elements are not repeated.
    void putRootBuds(set<Bud>& tray) const;
    void putBud(set<Bud>& tray) const;

    void putQuadraticForm(SymMat<VCData>& Q) const;
//    bool putQuadraticForm(SymMat<Double>& Q, multimap<Int4, VecDat3<Int4> >& qindex_hkl) const;

    void print(ostream&, const Double& minQ, const Double& maxQ) const;

    template <class Iterator>
    static void print(const string&, const Iterator&, const Iterator&);
};


//inline Int4 TreeLattice::SuperBasis() const
//{
//	if( HeadIsTail ) return 2;
//	if( m_root != NULL && m_root_on_left != NULL && m_root_on_right != NULL ) return 1;
//	return 0;
//}


inline const NodeB& TreeLattice::Root() const
{
	if( m_root == NULL )
	{
       	throw ZErrorMessage(ZErrorNullPointer, __FILE__, __LINE__, __FUNCTION__);
	}
	return *m_root;
}


inline void TreeLattice::setRoot(const Int4& K1, const Int4& K2)
{
	this->clear();
	m_root = new NodeB(K1, K2);
}


inline void TreeLattice::setRoot(const NodeB& lhs, const NodeB& rhs)
{
	this->clear();
	assert( lhs.Right() == rhs.Right() );

	m_root = new NodeB(lhs, rhs);
}

inline void TreeLattice::setRootEqualUpper()
{
	delete m_root_on_left;
	m_root_on_left = NULL; 
	delete m_root_on_right;
	m_root_on_right = NULL; 

	if( m_root->IsBud() ) HeadIsTail = false;
	else HeadIsTail = true;
	m_is_set_sort_criteria = false;
}

inline void TreeLattice::setRootOnLeftBranch(const NodeB& nodex)
{
		assert( m_root != NULL );

		delete m_root_on_left;
		m_root_on_left = NULL; 

		if( nodex.Left() == m_root->Right() )
		{
			if( m_root_on_right == NULL )
			{
				HeadIsTail = false;
				m_root_on_left = new NodeB(nodex);
				m_root_on_right = new NodeB(m_root_on_left->Right(), m_root->Left());
			}
			else if( nodex.Right() == m_root_on_right->Left() )
			{
				m_root_on_left = new NodeB(nodex);
			}
			else assert( false );
		}
		else if( nodex.Right() == m_root->Right() )
		{
			if( m_root_on_right == NULL )
			{
				HeadIsTail = false;
				m_root_on_left = new NodeB(nodex);
				m_root_on_left->swapBranch();
				m_root_on_right = new NodeB(m_root_on_left->Right(), m_root->Left());
			}
			else if( nodex.Left() == m_root_on_right->Left() )
			{
				m_root_on_left = new NodeB(nodex);
				m_root_on_left->swapBranch();
			}
			else assert( false );
		}
		else assert( false );

	m_is_set_sort_criteria = false;
}


inline void TreeLattice::setRootOnRightBranch(const NodeB& nodex)
{
	assert( m_root != NULL );

	delete m_root_on_right;
	m_root_on_right = NULL;
		
	if( nodex.Right() == m_root->Left() )
	{
		if( m_root_on_left == NULL )
		{
			HeadIsTail = false;
			m_root_on_right = new NodeB(nodex);
			m_root_on_left = new NodeB(m_root->Right(), m_root_on_right->Left());
		}
		else if( nodex.Left() == m_root_on_left->Right() )
		{
			m_root_on_right = new NodeB(nodex);
		}
		else assert( false );
	}
	else if( nodex.Left() == m_root->Left() )
	{
		if( m_root_on_left == NULL )
		{
			HeadIsTail = false;
			m_root_on_right = new NodeB(nodex);
			m_root_on_right->swapBranch();
			m_root_on_left = new NodeB(m_root->Right(), m_root_on_right->Left());
		}
		else if( nodex.Right() == m_root_on_left->Right() )
		{
			m_root_on_right = new NodeB(nodex);
			m_root_on_right->swapBranch();
		}
		else assert( false );
	}
	else assert( false );

	m_is_set_sort_criteria = false;
}


inline bool operator<(const TreeLattice& lhs, const TreeLattice& rhs)
{
	if( lhs.putCountOfQ() > rhs.putCountOfQ() ) return true;
	if( lhs.putCountOfQ() < rhs.putCountOfQ() ) return false;
	return lhs.putAreaSquare() < rhs.putAreaSquare();
}

inline void TreeLattice::swapBranch()
{
	if( m_root != NULL ) m_root->swapBranch();
	swap(m_root_on_left, m_root_on_right);

	if( m_root_on_left != NULL ) m_root_on_left->swapBranch();
	if( m_root_on_right != NULL ) m_root_on_right->swapBranch();
}



template <class Iterator>
void TreeLattice::print(const string& fname,
		const Iterator& it_begin, const Iterator& it_end)
{
	const vector<QData>& qdata = VCData::putPeakQData();
	const Double maxQ = qdata.rbegin()->q; // *max_element(Qdata.begin(), Qdata.end());
	const Double minQ = qdata.begin()->q; // *min_element(Qdata.begin(), Qdata.end());

	ofstream ofs(fname.c_str());

	ofs << "** MaxQ = " << maxQ << endl;
	ofs << "** MinQ = " << minQ << endl;
	ofs << "**" << endl; 
	
	// Output
    Int4 index = 1;
    for(Iterator it = it_begin; it!=it_end; it++)
    {
    	ofs << "** Tree_" << index++ << endl; 
    	it->print(ofs, minQ, maxQ);
        ofs << endl;
    }
}


#endif

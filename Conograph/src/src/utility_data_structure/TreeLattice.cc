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
#include <algorithm>
#include <cmath>
#include "FracMat.hh"
#include "TreeLattice.hh"

// RelationMatrix TreeLattice::m_rel_mat;

TreeLattice::TreeLattice()
{
	m_root = NULL;
	m_root_on_left = NULL;
	m_root_on_right = NULL;
	HeadIsTail = false;
	m_is_set_sort_criteria = false;
	m_count_Q = 0;
	m_detS = 0.0;
}


TreeLattice::TreeLattice(const TreeLattice& rhs)
{
	m_root = NULL;
	m_root_on_left = NULL;
	m_root_on_right = NULL;

	if( rhs.m_root != NULL ) m_root = new NodeB(*rhs.m_root);
	HeadIsTail = rhs.HeadIsTail;
	if( rhs.m_root_on_left != NULL )
	{
		m_root_on_left = new NodeB(*rhs.m_root_on_left);
	}
	if( rhs.m_root_on_right != NULL )
	{
		m_root_on_right = new NodeB(*rhs.m_root_on_right);
	}
	m_is_set_sort_criteria = rhs.m_is_set_sort_criteria;
	m_count_Q = rhs.m_count_Q;
	m_detS = rhs.m_detS;
}


TreeLattice::~TreeLattice()
{
	delete m_root;
	m_root = NULL;

	delete m_root_on_left;
	m_root_on_left = NULL;

	delete m_root_on_right;
	m_root_on_right = NULL;
}

TreeLattice& TreeLattice::operator=(const TreeLattice& rhs)
{
	if (this != &rhs)
	{
		clear();

		if( rhs.m_root != NULL ) m_root = new NodeB(*rhs.m_root);
		HeadIsTail = rhs.HeadIsTail;
		if( rhs.m_root_on_left != NULL )
		{
			m_root_on_left = new NodeB(*rhs.m_root_on_left);
		}
		if( rhs.m_root_on_right != NULL )
		{
			m_root_on_right = new NodeB(*rhs.m_root_on_right);
		}
		m_is_set_sort_criteria = rhs.m_is_set_sort_criteria;
		m_count_Q = rhs.m_count_Q;
		m_detS = rhs.m_detS;
	}
	return *this;
}


void TreeLattice::clear()
{
	delete m_root;
	m_root = NULL;

	delete m_root_on_left;
	m_root_on_left = NULL;

	delete m_root_on_right;
	m_root_on_right = NULL;

	HeadIsTail = false;
	m_is_set_sort_criteria = false;
	m_count_Q = 0;
	m_detS = 0.0;
}

void TreeLattice::setCountOfQ()
{
	if( m_root == NULL )
	{
		m_count_Q = 0;
	}
	else
	{
		set<Int4> index_tray;
		m_root->count(index_tray);
		if( m_root_on_left != NULL )
		{
			m_root_on_left->count(index_tray);
		}
		if( m_root_on_right != NULL )
		{
			m_root_on_right->count(index_tray);
		}

		m_count_Q = index_tray.size();
	}
}


void TreeLattice::setAreaSquare()
{
	set<Bud> budtray;
	putRootBuds(budtray);

	assert( !( budtray.empty() ) );

	m_detS = budtray.begin()->cross_product_312();
}


void TreeLattice::putRootBuds(set<Bud>& budtray) const
{
	budtray.clear();
	if( m_root == NULL ) return;
	
	if( HeadIsTail )
	{
		m_root->putRootBud(m_root->Upper(), budtray);
	}
	else if( m_root_on_left != NULL )
	{
		m_root->putRootBud(m_root_on_left->Right(), budtray);
		m_root_on_left->putRootBud(m_root->Left(), budtray);
		if( m_root_on_right != NULL ) m_root_on_right->putRootBud(m_root->Right(), budtray);
	}
	else if( m_root_on_right != NULL )
	{
		m_root->putRootBud(m_root_on_right->Left(), budtray);
		m_root_on_right->putRootBud(m_root->Right(), budtray);
	}
	else
	{
		m_root->putRootBud(-1, budtray);
	}
}


void TreeLattice::putBud(set<Bud>& budtray) const
{
	budtray.clear();
	if( m_root == NULL ) return;

	if( HeadIsTail )
	{
		m_root->putBud(m_root->Upper(), budtray);
	}
	else if( m_root_on_left != NULL )
	{
		m_root->putBud(m_root_on_left->Right(), budtray);
		m_root_on_left->putBud(m_root->Left(), budtray);
		if( m_root_on_right != NULL ) m_root_on_right->putBud(m_root->Right(), budtray);
	}
	else if( m_root_on_right != NULL )
	{
		m_root->putBud(m_root_on_right->Left(), budtray);
		m_root_on_right->putBud(m_root->Right(), budtray);
	}
	else
	{
		m_root->putBud(-1, budtray);
	}
}


//bool TreeLattice::putQuadraticForm(SymMat<Double>& Q, multimap<Int4, VecDat3<Int4> >& qindex_hkl) const
//{
//	static const VecDat3<Int4> hkl100(1,0,0);
//	static const VecDat3<Int4> hkl010(0,1,0);
//	static const VecDat3<Int4> hkl_1_10(-1,-1,0);
//
//	qindex_hkl.clear();
//	if( m_root == NULL ) return false;
//
//	if(m_root->Left() >= 0)
//	{
//		qindex_hkl.insert( multimap<Int4, VecDat3<Int4> >::value_type( m_root->Left(), hkl100 ) );
//	}
//	if(m_root->Right() >= 0)
//	{
//		qindex_hkl.insert( multimap<Int4, VecDat3<Int4> >::value_type( m_root->Right(), hkl010 ) );
//	}
//
//	if( m_root_on_left != NULL )
//	{
//		if(m_root_on_left->Right() >= 0)
//		{
//			qindex_hkl.insert( multimap<Int4, VecDat3<Int4> >::value_type( m_root_on_left->Right(), hkl_1_10 ) );
//		}
//
//		m_root->putQuadraticForm(hkl100, hkl010, qindex_hkl);
//		m_root_on_left->putQuadraticForm(hkl010, hkl_1_10, qindex_hkl);
//		if( m_root_on_right != NULL ) m_root_on_right->putQuadraticForm(hkl_1_10, hkl100, qindex_hkl);
//	}
//	else if( m_root_on_right != NULL )
//	{
//		if(m_root_on_right->Left() >= 0)
//		{
//			qindex_hkl.insert( multimap<Int4, VecDat3<Int4> >::value_type( m_root_on_right->Left(), hkl_1_10 ) );
//		}
//
//		m_root->putQuadraticForm(hkl100, hkl010, qindex_hkl);
//		m_root_on_right->putQuadraticForm(hkl_1_10, hkl100, qindex_hkl);
//	}
//	else
//	{
//		m_root->putQuadraticForm(hkl100, hkl010, qindex_hkl);
//	}
//
//	if( qindex_hkl.size() < 3 ) return false;
//
//	const vector<QData>& qdata = VCData::putPeakQData();
//	NRMat<Int4> icoef(3,3);
//	NRVec<Double> qvalue(3);
//
//	multimap<Int4, VecDat3<Int4> >::const_iterator it = qindex_hkl.begin();
//	icoef[0][0] = it->second[0]*it->second[0];
//	icoef[0][1] = it->second[1]*it->second[1];
//	icoef[0][2] = it->second[0]*it->second[1]*2;
//	qvalue[0] = qdata[(it++)->first].q;
//	icoef[1][0] = it->second[0]*it->second[0];
//	icoef[1][1] = it->second[1]*it->second[1];
//	icoef[1][2] = it->second[0]*it->second[1]*2;
//	qvalue[1] = qdata[(it++)->first].q;
//	icoef[2][0] = it->second[0]*it->second[0];
//	icoef[2][1] = it->second[1]*it->second[1];
//	icoef[2][2] = it->second[0]*it->second[1]*2;
//	qvalue[2] = qdata[(it++)->first].q;
//
//	const FracMat inv_mat = FInverse3( icoef );
//	assert( Q.size() == 2 );
//
//	Q(0,0) = ( inv_mat.mat[0][0]*qvalue[0] + inv_mat.mat[0][1]*qvalue[1] + inv_mat.mat[0][2]*qvalue[2] ) / inv_mat.denom;
//	Q(1,1) = ( inv_mat.mat[1][0]*qvalue[0] + inv_mat.mat[1][1]*qvalue[1] + inv_mat.mat[1][2]*qvalue[2] ) / inv_mat.denom;
//	Q(0,1) = ( inv_mat.mat[2][0]*qvalue[0] + inv_mat.mat[2][1]*qvalue[1] + inv_mat.mat[2][2]*qvalue[2] ) / inv_mat.denom;
//
//	return true;
//}


void TreeLattice::putQuadraticForm(SymMat<VCData>& Q) const
{
	assert(Q.size() == 2);
	if( m_root == NULL ) return;
	if( m_root->IsBud() ) return;

	set<Bud> budtray;
	this->putRootBuds(budtray);
	Q(0,0) = budtray.begin()->Q1();
	Q(1,1) = budtray.begin()->Q2();
	Q(0,1) = ( Q(0,0) + Q(1,1) - budtray.begin()->Q3() )/2;
}


void TreeLattice::print(ostream& os, const Double& minQ, const Double& maxQ) const
{
	if( m_root == NULL ) return;
	if( m_root->IsBud() ) return;
	
    os << "* Number of used peak positions: " << this->putCountOfQ() << endl;
    SymMat<VCData> Q(2);
    this->putQuadraticForm(Q);

    const Double det_Inv_Q = 1.0/( Q(0,0).Value()*Q(1,1).Value() - Q(0,1).Value()*Q(0,1).Value() );
    SymMat<Double> Inv_Q(2);
    Inv_Q(0,0) =  Q(1,1).Value()*det_Inv_Q;
    Inv_Q(1,1) =  Q(0,0).Value()*det_Inv_Q;
    Inv_Q(0,1) = -Q(0,1).Value()*det_Inv_Q;

    Double a = sqrt(Q(0,0).Value());
    Double b = sqrt(Q(1,1).Value());
    os << "* (a*, b*, \\gamma*): (" << a << ", " << b << ", " << acos(Q(0,1).Value()/(a*b))*180.0/M_PI << ")\n";

    a = sqrt(Inv_Q(0,0));
    b = sqrt(Inv_Q(1,1));
    os << "* (a, b, \\gamma): (" << a << ", " << b << ", " << acos(Inv_Q(0,1)/(a*b))*180.0/M_PI << ")\n";

    if( HeadIsTail )
	{
        os.width(15);
       	if( m_root->Upper() < 0 ) os << -1;
       	else os << m_root->Upper() + 1;	// HeadIsTail is true.

       	os.width(5);
    	os << " > ";
    	if( m_root->Upper() < 0 ) m_root->print(os, 20, -1, maxQ);
    	else m_root->print(os, 20, m_root->Upper(), maxQ);	// HeadIsTail is true.
        os << endl;
	}
	else if( m_root_on_left != NULL || m_root_on_right != NULL )
    {
        os.width(15);
    	if( m_root_on_left != NULL )
    	{
    		if( m_root_on_left->Right() < 0 ) os << -1;
        	else os << m_root_on_left->Right() + 1;

        	os.width(5);
        	os << " > ";
        	m_root->print(os, 20, m_root_on_left->Right(), maxQ);
    	}
    	else
    	{
    		if( m_root_on_right->Left() < 0 ) os << -1;
        	else os << m_root_on_right->Left() + 1;

        	os.width(5);
        	os << " > ";
        	m_root->print(os, 20, m_root_on_right->Left(), maxQ);
    	}
        os << endl;	

    	if( m_root_on_left != NULL )
    	{
            os.width(15);
        	if( m_root->Left() < 0 ) os << -1;
        	else os << m_root->Left() + 1;
        	os.width(5);
        	os << " > ";

        	m_root_on_left->print(os, 20, m_root->Left(), maxQ);
            os << endl;
    	}

    	if( m_root_on_right != NULL )
    	{
        	os.width(15);
        	if( m_root->Right() < 0 ) os << -1;
        	else os << m_root->Right() + 1;
        	os.width(5);
        	os << " > ";

        	m_root_on_right->print(os, 20, m_root->Right(), maxQ);
            os << endl;
    	}
    }
    else
    {
    	os << endl;
        os.width(15);
        if( m_root->Left() >= 0 && m_root->Right() >= 0 &&  m_root->Upper() >= 0 )
    	{
        	Double ans = ( VCData::putPeakPos(m_root->Left()).q + VCData::putPeakPos(m_root->Right()).q ) * 2.0
						- VCData::putPeakPos(m_root->Upper()).q;
        	if(ans < minQ) os << "<MinQ";
        	else os << ans;
    	}
        else os << -1;

        os.width(5);
    	os << " > "; 
    	if( m_root->Upper() < 0 ) m_root->print(os, 20, -1, maxQ);
    	else m_root->print(os, 20, m_root->Upper(), maxQ);	// HeadIsTail is true.
        os << endl;	
    }
}

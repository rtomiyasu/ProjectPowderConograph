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
#include "NodeB.hh"
#include "../zerror_type/error_out.hh"

NodeB::NodeB()
{
	m_left = -1;
	m_right = -1;
	
	m_left_branch = NULL;
	m_right_branch = NULL;
}


NodeB::NodeB(const Int4& i, const Int4& j)
{
	m_left = i;
	m_right = j;
	
	m_left_branch = NULL;
	m_right_branch = NULL;
}


NodeB::NodeB(const NodeB& lhs, const NodeB& rhs)
{
	m_left = lhs.Left();
	m_right = rhs.Left();
	assert( lhs.Right() == rhs.Right() );
	m_left_branch = new NodeB(lhs);
	m_right_branch = new NodeB(rhs);
}


NodeB::NodeB(const NodeB& rhs)
{
	m_left = rhs.m_left;
	m_right = rhs.m_right;
	
	m_left_branch = NULL;
	m_right_branch = NULL;
	if( rhs.m_left_branch != NULL ) m_left_branch = new NodeB(*rhs.m_left_branch);
	if( rhs.m_right_branch != NULL )
	{
		if( rhs.m_left_branch == rhs.m_right_branch ) m_right_branch = m_left_branch;
		else m_right_branch = new NodeB(*rhs.m_right_branch);
	}
}


NodeB::~NodeB()
{
	delete m_left_branch;
	m_left_branch = NULL;

	delete m_right_branch;
}


NodeB& NodeB::operator=(const NodeB& rhs)
{
	if (this != &rhs)
	{
		m_left = rhs.m_left;
		m_right = rhs.m_right;
		
		this->cutBranch();
		if( rhs.m_left_branch != NULL ) m_left_branch = new NodeB(*rhs.m_left_branch);
		if( rhs.m_right_branch != NULL )
		{
			if( rhs.m_left_branch == rhs.m_right_branch ) m_right_branch = m_left_branch;
			else m_right_branch = new NodeB(*rhs.m_right_branch);
		}
	}
	return *this;
}


void NodeB::setBranch(const Int4& index)
{
	cutBranch();
	
	m_left_branch = new NodeB(m_left, index);
	if( m_left == m_right ) m_right_branch = m_left_branch;
	else m_right_branch = new NodeB(m_right, index);
}


void NodeB::cutBranch()
{
	delete m_left_branch;
	m_left_branch = NULL;

	delete m_right_branch;
	m_right_branch = NULL;
}


void NodeB::count(std::set<Int4>& index_tray) const
{
	if( m_left == m_right )
	{
		if( m_left >= 0 ) index_tray.insert(m_left); 
		if( IsBud() ) return;
		
		m_left_branch->count(index_tray);
	}
	else
	{
		if( m_left >= 0 ) index_tray.insert(m_left); 
		if( m_right >= 0 ) index_tray.insert(m_right); 
		if( IsBud() ) return;

		m_left_branch->count(index_tray);
		m_right_branch->count(index_tray);
	}
}



void NodeB::putRootBud(const Int4& K3, std::set<Bud>& budtray) const
{
	Int4 count = 0;
	Int4 K4 = -1;
	if( !( this->IsBud() ) ) K4 = this->Upper();

	if( K3 < 0 ) count++;
	if( m_left < 0 ) count++;
	if( m_right < 0 ) count++;
	if( K4 < 0 ) count++;
	
	if( count < 2 )
	{
		Bud budex;
		budex.setIndex( m_left, m_right, K3, K4 );
		budtray.insert(budex);
	}
	else if( !( this->IsBud() ) )
	{
		m_left_branch->putRootBud(m_right, budtray);
		if( m_left != m_right ) m_right_branch->putRootBud(m_left, budtray);
	}
}


void NodeB::putBud(const Int4& K3, std::set<Bud>& budtray) const
{
	Int4 count = 0;
	Int4 K4 = -1;
	if( !( this->IsBud() ) ) K4 = this->Upper();

	if( K3 < 0 ) count++;
	if( m_left < 0 ) count++;
	if( m_right < 0 ) count++;
	if( K4 < 0 ) count++;

	if( count < 2 )
	{
		Bud budex;
		budex.setIndex( m_left, m_right, K3, K4 );
		budtray.insert(budex);
	}
	if( IsBud() ) return;

	m_left_branch->putBud(m_right, budtray);
	if( m_left != m_right ) m_right_branch->putBud(m_left, budtray);
}


//void NodeB::putQuadraticForm(const VecDat3<Int4>& hkl_left,
//						const VecDat3<Int4>& hkl_right,
//						multimap<Int4, VecDat3<Int4> >& qindex_hkl) const
//{
//	if( !( this->IsBud() ) )
//	{
//		VecDat3<Int4> diff = hkl_left - hkl_right;
//		if( this->Upper() >= 0)
//		{
//			qindex_hkl.insert( multimap<Int4, VecDat3<Int4> >::value_type( this->Upper(), diff ) );
//		}
//		m_left_branch->putQuadraticForm(hkl_left*(-1), diff, qindex_hkl);
//		if( m_left != m_right ) m_right_branch->putQuadraticForm(hkl_right, diff, qindex_hkl);
//	}
//}


void NodeB::print(ostream& os,  const Int4& num,
		const Int4& K3, const Double& maxQ) const
{
	os.width(10);
	os << this->chToString();

	bool flag_maxQ = false;
	bool set_Q4 = false;
	Double ans = -100.0;
	if( m_left_branch == NULL || m_right_branch == NULL )
	{
		if( m_left != -1 && m_right != -1 && K3 != -1 )
		{
			ans = ( VCData::putPeakPos(m_left).q
						+ VCData::putPeakPos(m_right).q ) * 2.0
						- VCData::putPeakPos(K3).q;
			if( ans > maxQ ) flag_maxQ = true;
			set_Q4 = true;
		}
	}
	
	os.width(5);
	os << " -L- "; 
	if( m_left_branch != NULL )
	{
		m_left_branch->print(os, num+15, m_right, maxQ);
	}
	else
	{
		if( flag_maxQ ) os << ">MaxQ";
		else if( set_Q4 ) os << ans;
	}
	os << endl;

	os.width(num+15);
	os << " -R- "; 
	if( m_right_branch != NULL )
	{
		m_right_branch->print(os, num+15, m_left, maxQ);
	}
	else
	{
		if( flag_maxQ ) os << ">MaxQ";
		else if( set_Q4 ) os << ans;
	}
	os << endl;
}

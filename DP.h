// DP.h: interface for the DP class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_DP_H__956559E8_00AB_4862_B6DB_0D9689308948__INCLUDED_)
#define AFX_DP_H__956559E8_00AB_4862_B6DB_0D9689308948__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

//#include	<stdio>
#include	<stdlib.h>
#include	<vector>
#include <stdio.h>
using namespace std;

class DP  
{
public:
	double			m_alpha;			// scale parameter

	// ss
	vector<int>		m_NumClassN;				// n(k): Number of descendant from each ancestor
	vector<vector<int> >	m_NumClassL;		// l(t,k)
	vector<vector<int> >	m_NumClassLA[2];	// la(t,k)
//la[0] is count of number of haplotype belonging to class k have 0 at that position t .	
//la[1] is count of number of haplotype belonging to class k have 1 at that position t .	
	vector<int>		m_pDataIndex;		// 

	int				m_NumClassU[3];	

public:
	int AddClass( int val, int numT );
	int DeleteClass();
	int DeleteSS();
	int LoadData( const char* filename );
	DP();
	virtual ~DP();

};

#endif // !defined(AFX_DP_H__956559E8_00AB_4862_B6DB_0D9689308948__INCLUDED_)

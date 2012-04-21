// DP.cpp: implementation of the DP class.
//
//////////////////////////////////////////////////////////////////////

#include "DP.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

DP::DP()
{
	m_alpha = 1; //0.4;
}

DP::~DP()
{

}

int DP::LoadData(const char *filename)
{
	int  t;
	FILE *fp = fopen( filename, "r" );

	char oneline[5000];
	while ( !feof(fp) )
	{
		fgets( oneline, 5000, fp );
		for (t=0; ;t++)
		{
			
		}

	}


	fclose(fp);
	return 1;	
}

int DP::AddClass(int val, int numT)
{
	m_NumClassN.push_back(val);	

	vector<int> zeros( numT, 0 );
	m_NumClassL.push_back( zeros );
	m_NumClassLA[0].push_back( zeros );		
	m_NumClassLA[1].push_back( zeros );
	
	return 1;
}

int DP::DeleteSS()
{
	m_NumClassN.clear();	
	m_NumClassL.clear();
	m_NumClassLA[0].clear(); 
	m_A.clear;		
	m_NumClassLA[1].clear();
	
	return 1;
}

int DP::DeleteClass()
{
	if ( m_NumClassN.size() > 0 )
	{
		m_NumClassN.pop_back();	

		m_NumClassL.pop_back();
		m_NumClassLA[0].pop_back();		
		m_NumClassLA[1].pop_back();
	}
	
	return 1;
}

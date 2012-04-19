// HapAssembler.h: interface for the HapAssembler class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_HAPASSEMBLER_H__E778E265_70BA_4A94_9F96_99879A8745D2__INCLUDED_)
#define AFX_HAPASSEMBLER_H__E778E265_70BA_4A94_9F96_99879A8745D2__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <vector>
#include <string>
#include "GenoHaploDB.h"
#include "HDP.h"

using namespace std;

class HapAssembler  
{
	int			m_nNumIteration2;
	string		m_Input_Dir;
	string		m_InputName_Header;
	int			m_nStart_t;
	int			m_nLast_t;
	
	vector<vector<unsigned char> >	m_Anc;
	vector<vector<int> >			m_CandiPair;
	vector<int>						m_NumClassN;
	int				*m_EqClass[2];		// C_i1, C_i2
	vector<float**>					m_LocalFreq;

public:
	HDP*		m_pHaplo;

	string				m_strOutdir;
	GenoHaploDB*		m_pDB;
	vector<haplo2_t>	m_arrData;
	vector<haplo2_t>	m_arrPairwisePred;
	vector<float>		m_BlockConfidence;

public:

	int LigateHaplotypes( vector<haplo2_t> &arrData, 
		vector<int> &T1, vector<int> &nOverlap, 
		vector<haplo2_t> &arrPred );
	
	int HierLigation( vector<haplo2_t> &arrData, vector<int> &T1, vector<int> &nOverlap, 
					  vector<haplo2_t> &arrPred );
	int PairwiseLigation( vector<haplo2_t> &arrData, 
		vector<int> &T1, vector<int> &nOverlap, 
					vector<haplo2_t> &arrPred );
	int LigateAdj(unsigned char **hl[], unsigned char **hr[], 
			int I, int T1, int T2, int nOverlap, 
			unsigned char **hout[], haplo2_t hest = 0, int offset = 0 );

	void	Initialize( GenoHaploDB* pDB, HDP* pHaplo ) 
				{ m_pDB = pDB; m_pHaplo = pHaplo; 	} ;

	int SaveResult( char const* ofilename );	
	int LoadAll( int tstart, int tend, int bLength, 
			int nOverlap, int numI, vector<int> &T1, vector<int> &vOverlap );
	
	int AddAnc( unsigned char* hl, unsigned char* hr, 
						 int T1, int T2, int nOverlap, 
						 vector<vector<unsigned char> > &Anc, int pref = 0);
	int Find( vector<unsigned char> h, 
		vector<vector<unsigned char> > &A,   int exceptK=-1  );
	int Findvec( unsigned char* h, 
		vector<vector<unsigned char> > &A,   int exceptK=-1  );
	float	GetEntropy( vector<int> &count );

	HapAssembler();
	HapAssembler(string &indir, string &header);
	virtual ~HapAssembler();

};

#endif // !defined(AFX_HAPASSEMBLER_H__E778E265_70BA_4A94_9F96_99879A8745D2__INCLUDED_)

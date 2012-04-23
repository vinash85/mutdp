// GenoHaploDB.h: interface for the GenoHaploDB class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_GENOHAPLODB_H__A16AF999_B19C_44BD_8FE9_674DA1977213__INCLUDED_)
#define AFX_GENOHAPLODB_H__A16AF999_B19C_44BD_8FE9_674DA1977213__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <vector>
#include "program.h"

using namespace std;

#define		NumAllele	2

typedef unsigned char*** haplo2_t;

class GenoHaploDB  
{
public:
	// data description
	int m_nGroups;
	int	m_numTotalI;
	int	m_numT;
	int	m_nMinT;
	int	m_nMaxT;

	int m_nTotHeteroSites ;
	int m_nTotHeteroInd ;
	
	unsigned char	**m_TrueHaplotypes[2];
	unsigned char	**m_RawGenotypes;
	
	unsigned char	**m_Genotypes[2];
	
	vector<int>		m_DataInGroup;
	vector<int>		m_numI;


	// program variables
	unsigned char	**m_Haplotypes[2];
	int				**m_Cum_Haplotypes[B][2];
	unsigned char	**m_Pred_Haplotypes[2];

	int				*m_EthnicGroup;		//

	int		**m_g_match;
	int		**m_g_miss1;
	int		**m_g_miss2;
	int		**m_h_count[2];
	bool	**m_FromTopLevel;

	bool	m_bNoisyObservation;

public:
	int		CopyHaplo(int bs, int be, haplo2_t *phh);
	float EstimatedTheta( vector<int> &m, vector<vector<int> > &l,vector<float> &theta  );
	int UpdateCumHaplotype( unsigned char ***h, 
		unsigned char ***h_old, int nstart, int nend, 
		int **cross );
	template <class T>
	int DeAlloc2DMemory(T ***ppArray, int nr, int nc);
	template <class T>
	int Alloc2DMemory(T ***ppArray, int nr, int nc);
	int AllocHaplo2( haplo2_t *hh, int numI, int bLength );
	int AllocMemory( int numI, int numT );
	int LoadData( vector<string> &filelist );
	GenoHaploDB();
	virtual ~GenoHaploDB();

};




#endif // !defined(AFX_GENOHAPLODB_H__A16AF999_B19C_44BD_8FE9_674DA1977213__INCLUDED_)

// GenoHaploDB.cpp: implementation of the GenoHaploDB class.
//
//////////////////////////////////////////////////////////////////////
#include <vector>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include "GenoHaploDB.h"
template <class T>
extern T	sum( vector<T> &invec );
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

GenoHaploDB::GenoHaploDB()
{
	m_nGroups = 0;
	m_numTotalI = 0;
	m_nTotHeteroSites = 0;
	m_nTotHeteroInd = 0;

	m_bNoisyObservation = 0;
	m_nMaxT = 10000;
}

GenoHaploDB::~GenoHaploDB()
{
	
}

// load groups of data 
int GenoHaploDB::LoadData( vector<string> &filelist )
{
	int i, j, t;
	int nfile = filelist.size();

	m_nGroups = nfile;
	m_numTotalI = 0;
	m_DataInGroup = new vector<int>[ m_nGroups ];

	int	numI, numT;
	for (j=0; j < nfile; j++)
	{		
		FILE *fp = fopen( filelist[j].c_str(), "r" );
		if ( fp == 0 )
		{
			printf( "Invalid file name!: %s\n", filelist[j].c_str() );
			return 0;
		}
		
		fscanf( fp, "%d %d", &numI, &numT );
		m_numI.push_back(numI);
		m_numTotalI += numI;
		fclose(fp);
	}

	AllocMemory( m_numTotalI, numT );

	int iIndex = 0;
	for (j=0; j < nfile; j++)
	{
		const char *filename = filelist[j].c_str();
		FILE *fp = fopen( filename, "r" );
		
		int	numI, numT;
		fscanf( fp, "%d %d", &numI, &numT );
		
		// load True haplotypes
		int tmp;
		for (i=0; i<numI; i++)
		{
			for (t=0; t<numT; t++)
			{
				fscanf( fp, "%d ", &tmp );
				m_TrueHaplotypes[0][i+iIndex][t] = tmp;
			}
			for (t=0; t<numT; t++)
			{
				fscanf( fp, "%d ", &tmp );
				m_TrueHaplotypes[1][i+iIndex][t] = tmp;
			}
		}
		if ( !m_bNoisyObservation )
		{
			// Draw Genotypes from True haplotypes
			for (i=0; i<numI; i++)
			{
				for (t=0; t<numT; t++)
				{
					if ( m_TrueHaplotypes[0][i+iIndex][t] == m_TrueHaplotypes[1][i+iIndex][t] )
						m_RawGenotypes[i+iIndex][t] =  m_TrueHaplotypes[0][i+iIndex][t];
					else
						m_RawGenotypes[i+iIndex][t] = 2;			
				}
				
				for (t=0; t<numT; t++)
				{
					if ( m_TrueHaplotypes[0][i+iIndex][t] <= m_TrueHaplotypes[0][i+iIndex][t] )
					{
						m_Genotypes[0][i+iIndex][t] =  m_TrueHaplotypes[0][i+iIndex][t]; //+ 1;
						m_Genotypes[1][i+iIndex][t] =  m_TrueHaplotypes[1][i+iIndex][t]; // + 1;
					}
					else
					{
						m_Genotypes[0][i+iIndex][t] =  m_TrueHaplotypes[1][i+iIndex][t]; // + 1;
						m_Genotypes[1][i+iIndex][t] =  m_TrueHaplotypes[0][i+iIndex][t]; // + 1;
					}
				}
			}
		}
		else
		{
			string gfile( filelist[j].c_str() );
			gfile = gfile.substr(0, gfile.length()-1);
			gfile += 'g';
			FILE *fp = fopen( gfile.c_str(), "r" );
			
			int	nI, nT;
			fscanf( fp, "%d %d", &nI, &nT );
			
			// load genotypes directly from files
			int tmp;
			for (i=0; i<numI; i++)
			{
				for (t=0; t<numT; t++)
				{
					fscanf( fp, "%d ", &tmp );
					m_RawGenotypes[i+iIndex][t] = tmp;
				}

				for (t=0; t<numT; t++)
				{
					if ( m_RawGenotypes[i+iIndex][t] == 0 )
					{
						m_Genotypes[0][i+iIndex][t] =  0; //+ 1;
						m_Genotypes[1][i+iIndex][t] =  0; // + 1;
					}
					else if ( m_RawGenotypes[i+iIndex][t] == 1 )
					{
						m_Genotypes[0][i+iIndex][t] =  1; // + 1;
						m_Genotypes[1][i+iIndex][t] =  1; // + 1;
					}
					else
					{
						m_Genotypes[0][i+iIndex][t] =  0; // + 1;
						m_Genotypes[1][i+iIndex][t] =  1; // + 1;
					}
				}
			}
			fclose( fp );
		}

		// Ethnic group label
		for (i = 0; i < numI; i++)
		{
			m_DataInGroup[j].push_back(i+iIndex);
			m_EthnicGroup[i+iIndex] = j;
		}
		iIndex += numI;

	}

	m_numT = MIN( numT, m_nMaxT );

	return 0;
}



int GenoHaploDB::AllocMemory(int numI, int numT)
{
	// haplotypes
	Alloc2DMemory( &m_Haplotypes[0], numI, numT );
	Alloc2DMemory( &m_Haplotypes[1], numI, numT );

	for (int bb=0; bb<B; bb++)
	{
		Alloc2DMemory( &m_Cum_Haplotypes[bb][0], numI, numT );
		Alloc2DMemory( &m_Cum_Haplotypes[bb][1], numI, numT );
	}

	// pred_H
	Alloc2DMemory( &m_Pred_Haplotypes[0], numI, numT );
	Alloc2DMemory( &m_Pred_Haplotypes[1], numI, numT );
	
	Alloc2DMemory( &m_TrueHaplotypes[0], numI, numT );
	Alloc2DMemory( &m_TrueHaplotypes[1], numI, numT );

	// genotypes
	Alloc2DMemory( &m_RawGenotypes, numI, numT );

	Alloc2DMemory( &m_Genotypes[0], numI, numT );
	Alloc2DMemory( &m_Genotypes[1], numI, numT );
	
	m_EthnicGroup = new int[numI];
	memset( m_EthnicGroup, 0, sizeof(int)*numI );

	m_DataInGroup = new vector<int>[m_nGroups];
	
	Alloc2DMemory( &m_FromTopLevel, m_numTotalI, 2);
	Alloc2DMemory( &m_h_count[0], m_numTotalI, numT);
	Alloc2DMemory( &m_h_count[1], m_numTotalI, numT);
	Alloc2DMemory( &m_g_match, m_numTotalI, numT);
	Alloc2DMemory( &m_g_miss1, m_numTotalI, numT);
	Alloc2DMemory( &m_g_miss2, m_numTotalI, numT);

	return 1;
}

template <class T>
int GenoHaploDB::Alloc2DMemory(T ***ppArray, int nr, int nc)
{
	int i;

	*ppArray = new T*[nr];

	for (i=0; i<nr; i++)
	{
		(*ppArray)[i] = new T[nc];
		memset( (*ppArray)[i], 0, sizeof(T)*nc );
	}

	return 1;
}

template <class T>
int GenoHaploDB::DeAlloc2DMemory(T ***ppArray, int nr, int nc)
{
	int i;

	for (i=0; i<nr; i++)
	{
		delete [] (*ppArray)[i];
		(*ppArray)[i] = NULL;
	}

	delete [] (*ppArray);
	(*ppArray) = NULL;

	return 1;
}

int GenoHaploDB::AllocHaplo2( haplo2_t *hh, int numI, int bLength )
{
	*hh = new unsigned char**[2];

	Alloc2DMemory( &( (*hh)[0]), numI, bLength );
	Alloc2DMemory( &( (*hh)[1]), numI, bLength );
	
	return 1;
}

int GenoHaploDB::UpdateCumHaplotype(
	unsigned char ***h, unsigned char ***h_old,
	int nstart, int nend, int **cross )
{
	int I = m_numTotalI;
	int ii, tt, bb;

	if ( h_old == 0 )
	{
		for ( ii=0; ii<I; ii++)
		{
			for ( tt = nstart; tt < nend ; tt++ )
			{
				bb = h[0][ii][tt];
				m_Cum_Haplotypes[bb][0][ii][tt] = 1;
				m_Cum_Haplotypes[1-bb][0][ii][tt] = 0;
				bb = h[1][ii][tt];
				m_Cum_Haplotypes[bb][1][ii][tt] = 1;
				m_Cum_Haplotypes[1-bb][1][ii][tt] = 0;
			}
		}
	}
	else
	{
		for ( ii=0; ii<I; ii++)
		{
			int mat1=0, mat2=0;
			for ( tt = nstart; tt < nend ; tt++ )
			{
				if ( h[0][ii][tt] == h_old[0][ii][tt-nstart] )
					mat1++;
				if ( h[1][ii][tt] == h_old[1][ii][tt-nstart] )
					mat1++;
				if ( h[0][ii][tt] == h_old[1][ii][tt-nstart] )
					mat2++;
				if ( h[1][ii][tt] == h_old[0][ii][tt-nstart] )
					mat2++;
			}
			
			if ( mat1 < mat2 )
				(*cross)[ii] += 1;

			if ( (*cross)[ii] % 2 == 0 )
			{
				for ( tt = nstart; tt < nend ; tt++ )
				{
					bb = h[0][ii][tt];
					m_Cum_Haplotypes[bb][0][ii][tt] += 1;
					bb = h[1][ii][tt];
					m_Cum_Haplotypes[bb][1][ii][tt] += 1;
				}
			}
			else
			{
				for ( tt = nstart; tt < nend ; tt++ )
				{
					bb = h[1][ii][tt];
					m_Cum_Haplotypes[bb][0][ii][tt] += 1;
					bb = h[0][ii][tt];
					m_Cum_Haplotypes[bb][1][ii][tt] += 1;
				}

			}
		}
		
	}

	return 1;
}

float GenoHaploDB::EstimatedTheta( vector<int> &m, vector<vector<int> > &l,
								vector<float> &theta )
{
	theta.clear();
	int K = m.size();
	float Tf = l[0].size();

	float sumtheta = 0;
	float theta_k;
	int invalid  = 0;
	for (int kk=0; kk < K; kk++)
	{
		vector<int> &lk = l[kk];
		float sumlk = (float)sum( lk );
		if ( m[kk] != 0 )
		{
			theta_k = 1 - sumlk / m[kk] /Tf;
			sumtheta += theta_k;
		}
		else
		{
			theta_k = -1;
			invalid++;
		}
		theta.push_back( theta_k );
	}
	float thetamean = sumtheta / (float)(K-invalid);
	
	return thetamean;
}

int GenoHaploDB::CopyHaplo(int bs, int be, haplo2_t *phh)
{
        int ee, ii;
        haplo2_t hh = *phh;

        for (ee=0; ee<2; ee++)
        {
                for (ii=0; ii< m_numTotalI; ii++)
                {
                        memcpy( hh[ee][ii],  &(m_Pred_Haplotypes[ee][ii][bs]),
                                (be-bs)*sizeof(unsigned char) );
                }
        }

        return 1;
}

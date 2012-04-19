// HapAssembler.cpp: implementation of the HapAssembler class.
//
//////////////////////////////////////////////////////////////////////

#include <algorithm> 
#include <string> 
#include <vector> 
#include <memory.h> 
#include <math.h>
#include <assert.h>

#include "HapAssembler.h"
#include "GenoHaploDB.h"
#include "HDP.h"
//#include "HapDP.h"

extern int normalize( vector<double> &vec );
template <class T>
extern int sample_discrete( vector<T> &pred );
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////


HapAssembler::HapAssembler()
{
	m_Input_Dir = string("./output/");
	m_InputName_Header = "Hapmap2_pred_";

	m_EqClass[0] = 0;
	m_EqClass[1] = 0;
	m_nNumIteration2 = 3000;

	m_strOutdir = string("./output/");
	m_pHaplo = 0;
	
}

HapAssembler::HapAssembler(string &indir, string &header)
{
	m_Input_Dir = indir;
	m_InputName_Header = header;
		
	m_EqClass[0] = 0;
	m_EqClass[1] = 0;
	m_nNumIteration2 = 3000;
	
	m_strOutdir = string("./output/");
	m_pHaplo = 0;
}

HapAssembler::~HapAssembler()
{

}

////////////////////////////////////////////////////////////
// ligate two haplotype blocks with I individual
// No error model in inheritance
// Retype if the diversity is not small enough
// -> improvementpoint: when retyping, make use of the information

// of the 'determined individuals'
int HapAssembler::LigateAdj(unsigned char **hl[], unsigned char **hr[], 
			int I, int T1, int T2, int nOverlap, 
			unsigned char **hout[], haplo2_t hest, int offset )
{
	int nT = T1 + T2 - nOverlap;
	int i, ii, tt, j, ee;

	m_Anc.clear();
	m_CandiPair.clear();
	m_NumClassN.clear();

	if ( m_EqClass[0] ) delete [] m_EqClass[0];
	if ( m_EqClass[1] ) delete [] m_EqClass[1];

	//equivalence class for father and mother
	//m_Eqclass give current haplotype belongs to which ancestor
	// i.e.  Ancestor(m_EqClass) = current haplotype
	m_EqClass[0] = new int[I];
	m_EqClass[1] = new int[I];


	//////////////////////////////////////////////////
	bool	*bDone = new bool[I];
	memset( bDone, false, I*sizeof(bool) );

	int diff1 = 0, diff2 = 0;
    	int kk[4];
	vector<double> Prob;
	vector<int> pair;
	int numDone = 0;
	/////////////////////////////////////////////
	// initialization
	for ( i=0; i < I; i++ )
	{
		bool bhetero = false;
		// see if homogeneous 
		int nhet = 0;
		for ( tt=0; tt< nOverlap; tt++ )
		{
			if ( hr[0][i][tt] != hr[1][i][tt] )
			{
				bhetero = true;				
				nhet++;
			}	
		}

		diff1 = 0; diff2 = 0;
		for ( tt=0; tt< nOverlap; tt++ )
		{
			if ( hl[0][i][tt+(T1-nOverlap)] != hr[0][i][tt] )
				diff1 ++;
			if ( hl[1][i][tt+(T1-nOverlap)] != hr[1][i][tt] )
				diff1 ++;
			if ( hl[0][i][tt+(T1-nOverlap)] != hr[1][i][tt] )
				diff2 ++;
			if ( hl[1][i][tt+(T1-nOverlap)] != hr[0][i][tt] )
				diff2 ++;
		}

		pair.clear();	
		// have heterogygous sites in the overlapping region
		if ( bhetero && ( diff1==0 || diff2==0 ) && ( diff1 != diff2)  ) 
		{	
			// compute the dissimilarity on overlapping regions
			bDone[i] = true;
			if (  diff1 == 0  ) 
			{
				kk[0] = AddAnc( hl[0][i], hr[0][i], T1, T2, nOverlap, m_Anc );
				kk[1] = AddAnc( hl[1][i], hr[1][i], T1, T2, nOverlap, m_Anc );
			}
			else
			{
				kk[0] = AddAnc( hl[0][i], hr[1][i], T1, T2, nOverlap, m_Anc );
				kk[1] = AddAnc( hl[1][i], hr[0][i], T1, T2, nOverlap, m_Anc );
			}
			int added = MAX( kk[0], kk[1] ) - m_NumClassN.size() + 1;
			while (added > 0 )
			{
				m_NumClassN.push_back(0);
				added --;
			}
			for ( ee = 0; ee < 2; ee++ )
			{
				m_EqClass[ee][i] = kk[ee];
				m_NumClassN[kk[ee]] += 1;
				pair.push_back( kk[ee] );
			}
			m_CandiPair.push_back( pair );
		}
		else
		{
			float randsel = MAX((float)rand()/RAND_MAX, 0.00001);
			for ( int rt=0;  rt <= nOverlap  && bDone[i]==false ; rt++ )
			{
				kk[0] = AddAnc( hl[0][i], hr[0][i], T1, T2, nOverlap, m_Anc, rt );
				kk[1] = AddAnc( hl[1][i], hr[1][i], T1, T2, nOverlap, m_Anc, rt );
				kk[2] = AddAnc( hl[0][i], hr[1][i], T1, T2, nOverlap, m_Anc, rt );
				kk[3] = AddAnc( hl[1][i], hr[0][i], T1, T2, nOverlap, m_Anc, rt );

				int added = MAX( MAX(kk[0], kk[1]), MAX(kk[2],kk[3]) ) - m_NumClassN.size() + 1;
				while (  added >  0 )
				{
					m_NumClassN.push_back(0);
					added --;
				}
				bDone[i] = false;
			
				for ( j = 0; j < 4; j+=2 )
				{
					bool bdup = false;
					for ( int ll=0; ll < pair.size(); ll+=2 )
					{
						if ( pair[ll] == kk[j] && pair[ll+1] == kk[j+1] )
						{
							bdup = true;
							break;
						}
					}
					if ( !bdup )
					{
						pair.push_back( kk[j] );
						pair.push_back( kk[j+1] );
					}
				}
			}
			if ( !bDone[i] )
			{
				m_CandiPair.push_back( pair );
				/// initialize with input estimation
				if ( hest )
				{
					m_EqClass[0][i] = Findvec( hest[0][i], m_Anc );
					m_EqClass[1][i] = Findvec( hest[1][i], m_Anc );
				}
				else
				{
					m_EqClass[0][i] = pair[0];
					m_EqClass[1][i] = pair[1];
				}
			}

			m_NumClassN[ m_EqClass[0][i] ] += 1;
			m_NumClassN[ m_EqClass[1][i] ] += 1;		
		}	// end of bhetero
		
		if ( bDone[i] ) numDone++;
	}


	int KK0 = m_NumClassN.size();
	int KK = 0;

	vector<vector<int> >	m_CumC1(I, vector<int>(m_Anc.size(), 0) );
	vector<vector<int> >	m_CumC2(I, vector<int>(m_Anc.size(), 0) );
	int nIter = m_nNumIteration2;
	int nBurnIn = nIter/2;

	for ( int iter = 0; iter < nIter; iter++ )
	{
		vector<int> iivec( I );
		for ( ii = 0; ii < I; ii++ )	iivec[ii] = ii;
		random_shuffle( iivec.begin(), iivec.end() );
	
		KK = 0;
		for (int k=0; k<KK0; k++)
		{
			if ( m_NumClassN[k] > 0 )
				KK++;
		}

		for ( ii=0; ii < I; ii++ )
		{
			i = iivec[ii];
//			if ( !bDone[i] ) 
			{
				int k1 = m_EqClass[0][i];
				int k2 = m_EqClass[1][i];

				m_NumClassN[k1]--;
				m_NumClassN[k2]--;

				vector<int> &candi = m_CandiPair[i];
				
				Prob.clear();

				for (j = 0; j < candi.size(); j+=2 )
				{
					Prob.push_back( ( m_NumClassN[ candi[j] ] + 1  ) *						 
						(m_NumClassN[ candi[j+1] ] + 1 ) );
				}
				normalize( Prob );
				int pk = 2*sample_discrete( Prob );

				m_EqClass[0][i] = candi[ pk ];
				m_EqClass[1][i] = candi[ pk+1 ];


				if ( m_EqClass[0][i] > m_EqClass[1][i] )
				{
					int tmp = m_EqClass[0][i];
					m_EqClass[0][i] = m_EqClass[1][i];
					m_EqClass[1][i] = tmp;
				}
				
				m_NumClassN[ m_EqClass[0][i] ] ++;
				m_NumClassN[ m_EqClass[1][i] ] ++;	

				if ( iter > nBurnIn )
				{
					m_CumC1[i][ m_EqClass[0][i] ] ++;	
					m_CumC2[i][ m_EqClass[1][i] ] ++;		
				}
			}		
		}
	}

	///////////////////////////////////////////////////////
	//////////////		posterior mean		///////////////
	///////////////////////////////////////////////////////
	for ( i=0; i < I; i++ )
	{
		int k1, k2;
		if (0) //( bDone[i] ) 
		{
			k1 = m_EqClass[0][i];
			k2 = m_EqClass[1][i];
			copy( m_Anc[ k1 ].begin(), m_Anc[ k1 ].end(), hout[0][i] );
			copy( m_Anc[ k2 ].begin(), m_Anc[ k2 ].end(), hout[1][i] );
		}
		else
		{
		   vector<int>::iterator max1 = 
			max_element( m_CumC1[i].begin(), m_CumC1[i].end());
		   vector<int>::iterator max2 = 
			max_element( m_CumC2[i].begin(), m_CumC2[i].end());

			k1 = (int) (max1 - m_CumC1[i].begin());
			k2 = (int) (max2 - m_CumC2[i].begin());
			copy( m_Anc[ k1 ].begin(), m_Anc[ k1 ].end(), hout[0][i] );
			copy( m_Anc[ k2 ].begin(), m_Anc[ k2 ].end(), hout[1][i] );
		}
	}

	haplo2_t h0;
	m_pDB->AllocHaplo2( &h0, I, nT );
	for (ee=0; ee<2; ee++)
	for (ii=0; ii<I; ii++)
		memcpy( h0[ee][ii], hout[ee][ii], sizeof(unsigned char)*nT );
	
	m_pHaplo->SetRange( 0, nT );
	m_pHaplo->GetPredFreq( h0 );
	float ent = GetEntropy( m_pHaplo->GetSumClassN() );

	if ( (ent > 2.5) && nT < 22 ) 
	{
		m_pHaplo->InferHaplotypes( h0, m_pDB->m_EthnicGroup, I, nT, 
			offset );
		
		m_pHaplo->Sample_Pred();
		for (int ee=0; ee<2; ee++)
		for (int ii=0; ii<I; ii++)
		for (int tt=0; tt<nT; tt++)
			hout[ee][ii][tt] = m_pDB->m_Pred_Haplotypes[ee][ii][tt+offset];

	}

	if ( bDone ) delete [] bDone;

	return 1;
}

int HapAssembler::PairwiseLigation( vector<haplo2_t> &arrData, 
					vector<int> &T1, vector<int> &nOverlap,
					vector<haplo2_t> &arrPred )
{
	int I = m_pDB->m_numTotalI;
	int numT = m_pDB->m_numT;
	unsigned char ***hRight;
	unsigned char ***hLeft;

	////// clear output array /////
	arrPred.clear();

	int numBlocks = arrData.size();

	int Tmax = 0;
	for (int bb=0; bb<numBlocks; bb++)
	{
		if ( T1[bb] > Tmax )
			Tmax = T1[bb];
	}
	m_pDB->AllocHaplo2( &hLeft, I, Tmax );
	m_LocalFreq.clear();

	int Tleft = T1[0];
	///		copy data to local place 
	for (int i=0; i<I; i++)
	{
		memcpy( hLeft[0][i], arrData[0][0][i], T1[0] );
		memcpy( hLeft[1][i], arrData[0][1][i], T1[0] );
	}

	if ( numBlocks > 1 )
	{
		int offset = 0;
		vector<int> T1new;
		vector<int> nOvernew;
		for ( int nb = 1 ; nb < numBlocks ; nb++ )
		{
			int bs = offset + T1[nb-1] - nOverlap[nb-1] ;
			int be = MIN( bs + T1[nb],  numT );

			hRight = arrData[nb];

			T1new.push_back( T1[nb-1] + T1[nb] - nOverlap[nb-1] );
			if ( nb == numBlocks-1 )
				nOvernew.push_back( 0 );
			else
				nOvernew.push_back( T1[nb] );
				
			haplo2_t hOut;
			m_pDB->AllocHaplo2( &hOut, I, T1new[nb-1] );
		
			///////////////////////////////////////
			///////// * do ligation ///////////////////
			LigateAdj( hLeft, hRight, I, T1[nb-1], T1[nb], nOverlap[nb-1], hOut, 0, offset );
			arrPred.push_back( hOut );
			
			//// copy hRight to hLeft for next iteration
			for (int i=0; i<I; i++)
			{
				memcpy( hLeft[0][i], hRight[0][i], T1[nb] );
				memcpy( hLeft[1][i], hRight[1][i], T1[nb] );
			}

			////////////////////////////////////
			///		m_pHaplo: just for some utility function. 
			/////////////////////////////////

			m_pHaplo->SetRange(  0, T1new[nb-1]  );	// adjust range to PredFreq
			m_pHaplo->GetPredFreq( hOut );			// adjust ss
			////////////////////////////////////	
			
			// return to the original range and compute error 
			m_pHaplo->SetRange( bs - (T1[nb-1] - nOverlap[nb-1]), be );

			fflush( stdin );

			if ( be == numT )
				break;
			
			offset += T1[nb-1] - nOverlap[nb-1];
		}

		T1.clear();
		T1 = T1new;
		nOverlap.clear();
		nOverlap = nOvernew;
	}
	else
	{
		arrPred.clear();
		arrPred.push_back( hLeft );
	}

	return 1;
}

int HapAssembler::LigateHaplotypes(vector<haplo2_t> &arrData, 
				vector<int> &T1, vector<int> &nOverlap, 
				vector<haplo2_t> &arrPred )
{

	//////////////////////////////////////////////
	/// 1. firstlevel pairwise ligation 
	PairwiseLigation( arrData, T1, nOverlap, arrPred );
	arrData.clear();

	////////// try once step again (optional)/////////////	
	arrData = arrPred;
	arrPred.clear();
	PairwiseLigation( arrData, T1, nOverlap, arrPred );
	arrData.clear();
	///////////////////////////////////

	int* popLabel = m_pDB->m_EthnicGroup;
	
	arrData = arrPred;
	
	arrPred.clear();
	
	//////////////////////////////////////////////
	///	2. hierarchical ligation 
	HierLigation( arrData, T1, nOverlap, arrPred );

	//////////////////////////////////////////////
	/// copy result
	int I = m_pDB->m_numTotalI;
	int numT = m_pDB->m_numT;
	haplo2_t &hh = arrPred[0];
	for (int i=0; i<I; i++)
	{
		memcpy( m_pDB->m_Pred_Haplotypes[0][i], (unsigned char*)hh[0][i], 
			sizeof(unsigned char)*numT );
		memcpy( m_pDB->m_Pred_Haplotypes[1][i], (unsigned char*)hh[1][i], 
			sizeof(unsigned char)*numT );
	}

	return 1;
}

int HapAssembler::AddAnc( unsigned char* hl, unsigned char* hr, 
						 int T1, int T2, int nOverlap, 
						 vector<vector<unsigned char> > &Anc, int roffset )
{
	// add h_l + h_r to the pool
	int nT = T1 + T2 - nOverlap;
	int k;

	vector<unsigned char> hnew( nT );
	unsigned char* pp = new unsigned char[nT];
	unsigned char* pTmp;
	if ( roffset == 0 )
	{
		pTmp = copy( hl, (hl+T1), pp );
		copy( &(hr[nOverlap]), &(hr[nOverlap]) + T2-nOverlap, pTmp );
	}
	else	
	{
		pTmp = copy( hl, hl+(T1-roffset), pp );
		copy( hr+(nOverlap-roffset), hr + T2, pTmp );
	}

	for (int i=0; i<nT; i++)
	{
		hnew[i] = pp[i];
	}
	k = Find( hnew, Anc );
	if ( k < 0 )
	{
		Anc.push_back( hnew );
		k = Anc.size() - 1;
	}
	if ( pp ) delete [] pp;

	return k;
}

int HapAssembler::Find( vector<unsigned char> h, 
	vector<vector<unsigned char> > &A,   int exceptK  )
{
	if ( A.size() == 0 )
		return -1;

	int		kk, tt;
	int		T = A[0].size();
	int		K = A.size();
	int		found = 0;
	
	for (kk=0; kk<K; kk++)
	{
		if ( kk != exceptK )
		{
			for (tt=0; tt<T; tt++)
			{
				if ( h[tt] != A[kk][tt] )
					break;
			}
			if ( tt == T )
				break;
		}
	}

	if ( kk < K )
		return kk;
	else
		return -1;
}

int HapAssembler::Findvec( unsigned char* h, 
	vector<vector<unsigned char> > &A,   int exceptK  )
{
	if ( A.size() == 0 )
		return -1;

	int		kk, tt;
	int		T = A[0].size();
	int		K = A.size();
	int		found = 0;
	
	for (kk=0; kk<K; kk++)
	{
		if ( kk != exceptK )
		{
			for (tt=0; tt<T; tt++)
			{
				if ( h[tt] != A[kk][tt] )
					break;
			}
			if ( tt == T )
				break;
		}
	}

	if ( kk < K )
		return kk;
	else
		return -1;
}

int HapAssembler::HierLigation( vector<haplo2_t> &arrData, 
				vector<int> &T1, vector<int> &nOverlap,
				vector<haplo2_t> &arrPred )
{
	int I = m_pDB->m_numTotalI;
	int numT = m_pDB->m_numT;
	unsigned char ***hRight;
	unsigned char ***hLeft;
	unsigned char ***htmp;
	
	m_pDB->AllocHaplo2( &htmp, I, numT );
	
	vector<int> T1new;
	vector<int> nOvernew;

	int numBlocks = arrData.size();
	
	m_LocalFreq.clear();

	int Tleft = T1[0];

	int bs = 0, be, be2;
	while ( numBlocks > 1 )
	{
		arrPred.clear();
		bs = 0;
		for ( int nb = 0 ; nb < numBlocks ; nb+=2 )
		{
			be = MIN( bs + T1[nb], numT );
			be2 = MIN( be + T1[nb+1] - nOverlap[nb], numT );

			if ( nb == numBlocks - 1 )
			{
				arrPred.push_back( arrData[nb] );
				T1new.push_back( be2 - bs );
				nOvernew.push_back( 0 );
				break;
			}
			
			T1new.push_back( be2 - bs );
			if ( nb < numBlocks - 1 )
				nOvernew.push_back( nOverlap[nb+1] );
			else
				nOvernew.push_back( 0 );

			hLeft = arrData[nb];
			hRight = arrData[nb+1];
		
			haplo2_t hOut;
			m_pDB->AllocHaplo2( &hOut, I, T1[nb] + T1[nb+1] - nOverlap[nb] );
			
			LigateAdj( hLeft, hRight, I, T1[nb], T1[nb+1], nOverlap[nb], hOut, 0, bs );
			arrPred.push_back( hOut );
			
			m_pHaplo->SetRange( bs, be2 );

			fflush( stdin );

			if ( be == numT )
				break;

			bs += (T1[nb] + T1[nb+1] - nOverlap[nb] - nOverlap[nb+1] ) ;
		}

		arrData.clear();
		arrData = arrPred;
		numBlocks = arrData.size();
		T1.clear();
		T1 = T1new;
		nOverlap.clear();
		nOverlap = nOvernew;
		
		T1new.clear();
		nOvernew.clear();
	}

	arrPred.clear();
	arrPred.push_back( arrData[0] );
	
	return 1;
}

int HapAssembler::SaveResult( char const* ofilename )
{
	int ii, tt, ee;
	FILE *fp = fopen( ofilename, "w" );

	// Save H
	int I = m_pDB->m_numTotalI;
	unsigned char **pred_h[2];
	pred_h[0] = m_pDB->m_Pred_Haplotypes[0];
	pred_h[1] = m_pDB->m_Pred_Haplotypes[1];

	for (ii=0; ii < I; ii++)
	{
		fprintf( fp, "%d, G%d\n", ii, m_pDB->m_EthnicGroup[ii] );
		for (ee=0; ee<2; ee++)
		{
			for (tt = m_pHaplo->m_nBlockStart; tt< m_pHaplo->m_nBlockEnd; tt++)
				fprintf( fp, "%d", pred_h[ee][ii][tt] );

			fprintf( fp, "\n" );
		}
	}

	fclose(fp);
	return 1;
}


int HapAssembler::LoadAll( int tstart, int tend, int bLength, 
		int nOverlap, int numI, 
		vector<int> &T1, vector<int> &vOverlap )
{
	int b;
	vector<haplo2_t>& arrData = m_arrData;
	int nBlock = (int)floor( 
		(float)(tend - tstart - nOverlap)/( bLength - nOverlap) ) + 1; 

	int go = bLength - nOverlap;
	string filename( m_Input_Dir );
	filename += m_InputName_Header;
	printf("filename=%s\n", filename.c_str());
	string bfilename ;
	char tmp[100];
	int ts, te;
	
	m_pHaplo->SetRange( 0, tend - tstart );
	arrData.clear();
	m_BlockConfidence.clear();

	T1.clear();
	vOverlap.clear();
	float maxee = -100000;
	for (b = 0; b < nBlock; b++ )
	{
		haplo2_t	hh;
		
		GenoHaploDB DB;
		DB.AllocHaplo2( &hh, numI, bLength );

		bfilename = filename; 
		ts = go * b;
		te = ts + bLength;
		if ( b > 0 )
		{
			if ( te <= m_pDB->m_numT )
			{
				vOverlap.push_back( nOverlap );
			}
			else  // ( te > m_pDB->m_numT )
			{
				vOverlap.push_back( nOverlap + te - m_pDB->m_numT );
				te = m_pDB->m_numT;
				ts = te - bLength;

			}
		}
		T1.push_back( bLength );

		sprintf( tmp, "T%d_%d.txt", ts, te-1 );
		bfilename += string( tmp );
		m_pHaplo->SetRange( 0, bLength );

		m_pHaplo->LoadOutput( bfilename.c_str(), 0, bLength, hh );

		m_pHaplo->SetRange( 0, bLength );
		
		float ee = GetEntropy( m_pHaplo->GetSumClassN() );
		m_BlockConfidence.push_back( ee );
		if ( ee > maxee )
			maxee = ee;
		arrData.push_back( hh );

	}

	vOverlap.push_back( 0 );
	
	return 1;
}

float HapAssembler::GetEntropy( vector<int> &count  )
{
	float val = 0;
	int K = count.size();
	float total = 0;
	int I2 = m_pDB->m_numTotalI * 2;


	for ( int kk = 0; kk < K; kk++)
	{
		float pi = (float)count[kk]/I2;
		if ( pi > 0 )
			val -= pi * log(pi);
	}

	return val;
}

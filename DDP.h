// HDP.h: interface for the HDP class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_HDP_H__B8C60B76_E5B0_4C89_B1A4_86D1ABBF602C__INCLUDED_)
#define AFX_HDP_H__B8C60B76_E5B0_4C89_B1A4_86D1ABBF602C__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <vector>
#include <memory.h>
#include "program.h"
#include "GenoHaploDB.h"

using namespace std;

class MUTDP
{

	GenoHaploDB	*m_pData;		// pointer to data 

	double m_gamma;				// bottom level scale parameter


	// ss 
	/*vector<int>		m_NumClassN;*/
	/*vector<bool>	m_Remove_Class;*/

	// for convergence check
	vector<int>		m_traceDiff;

	int		m_bAmbiguous;

	////////////// I/O /////////////
	FILE	*m_pFPTheta;
	string	m_strOutdir;
	string	m_strFilename;
	//******ADDED FROM DP ******// 
	//TODO scale parameter

	double			m_alpha_a;			// top level parameter

	// ss top level
	vector<int>		m_NumClassQ;				// q(k)=sum_j I(c(j) =k) : Number of descendant B from each ancestor A
	//keeps track of difference and similarities between A and B.
	vector<vector<int> >	m_NumClassM;		// m(k,t) = \sum_j I(b(j,t) =a(k,t) ) I( c(j) ==k)
	// As we are not sampling A there is no need of numClassMA

	// ss bottom level
	double			m_alpha_b;			// bottom level scale parameter
	vector<int>		m_NumClassP;				// p(j)= sum_i I(d(j) =j): Number of descendant H from each ancestor B 
	//keeps track of difference and similarities between A and B.
	vector<vector<int> >	m_NumClassL;		// l(j,t) = \sum_i I(h(i,t) =b(j,t) ) I( d(i) ==j)

	//la[0] is count of number of haplotype belonging to class k have 0 at that position t .
	//la[1] is count of number of haplotype belonging to class k have 1 at that position t .
	vector<vector<int> >	m_NumClassLA[2];	// la(j,t) =\sum_i  I(h_it = 0) I(d_i =j) ---- required for sampling of B.

	//TODO
	//int				m_NumClassU[3];


public:
	// random variables to sample 
	/*vector<vector<unsigned char> >	m_A;	// ancestral haplotypes*/
	/*int				**m_EqClass;		// individual haplotypes */

	/*vector<vector<unsigned char> >	m_Buffer_A;*/

	// local block definition 
	int		m_nBlockStart;		// block-starting SNP index 
	int		m_nBlockEnd;		// block-ending SNP index 
	int		m_nBlockLength;		// block length 

	// iteration counter
	int		m_nthIter;
	int		m_nBurninIteration;
	int		m_nCumIteration;
	int		m_nThining;

	// flags 
	int		m_doConparam;
	int		m_bCheckConvg;


	///avinash

	vector<vector<unsigned char> >	m_A;	// ancestral haplotypes
	int		**m_EqClass;		// individual haplotypes

	vector<vector<unsigned char> >	m_Buffer_A;
	vector<bool>	m_Remove_Class;

	//******ADDED FROM DP ENDS******//



	///////////////////////////////////////////////////////
	///		wrapper function for block-wise hap inference
	int		InferHaplotypes( haplo2_t h0, int* popLabel, int I, int T, int offset, bool* bDone = 0);
	int		InferHaplotypes( GenoHaploDB* pDB, int nstart, int nend );
	int 	InferHaplotypesAdaptive( GenoHaploDB* pDB, int nstart, int nend );	

	///////////////////////////////////////////////////////
	///		internal functions to do hap inf iterations 
	int		Init( GenoHaploDB *pDB, int nstart, int nend );
	int		Initialize(haplo2_t h0, int I, int T, 
			int offset = 0, bool bCopyShiftedRand = 0 );
	int		Initialize( GenoHaploDB *pDB, int nstart, int nend );
	int		Iterate_det_Gibbs_Met( int numIter, bool* bDone = 0 );
	int		Iterate_cum_Gibbs_Met( int numIter, bool* bDone = 0 );
	// do prediction
	int		Sample_Pred();

	/////	internal functions for hap inf 
	vector<unsigned char> Sample_A( int cc, unsigned char *h, bool new_class );
	int		Sample_H( unsigned char *h, unsigned char *h1, 
			unsigned char *g0, unsigned char *g1,
			vector<unsigned char> &ak, 
			int	*g_match_i, int *g_miss1, int *g_miss2,
			int nk, vector<int> &lk, vector<int> &lak0, vector<int> &lak1, 
			int *h_count, 
			int *u, int I );	
	int		Sample_EqClass(  unsigned char *h,  vector<vector<int> > &l, 
			double alpha0);//, bool *FromTopLevel );
	int		Sample_EqClass_Init(  unsigned char *h, vector<vector<int> > &l, 
			double alpha0);//, bool *FromTopLevel );
	int		Sample_Conparam( int numiter_a, int numiter_b  ) ;
	bool	TestAcceptance(int old_c, int new_c, unsigned char *h, 
			vector<unsigned char> &old_a, vector<unsigned char> &temp_a, 
			vector<int> &n, vector<vector<int> > &l );

	int		EstimateTheta( int iter  );

	// util functions 
	int		MergeList();
	int		CondenseList( int bBackupA = 0);
	int		ResetCumH();
	int		CalNumClassU();
	int		BackupOldSS( int old_c, vector<bool> &remove_class );
	int		BackupUpdateSS( int new_c,vector<bool> &remove_class, bool classAdded );
	int		DeleteSS( int ii, int ee, int cc );
	int		RollBack( int ii, int ee,  int old_c, int new_c );
	int		AddClass( int initvalue, int numT );

	int				clearSS();
	void			SetRange( int bs, int be ) 
	{ m_nBlockStart = bs; m_nBlockEnd = be; m_nBlockLength = be - bs; }
	double			GetGamma()	{ return	m_gamma; }
	void			SetGamma(double gg)	{ m_gamma = gg; }
	GenoHaploDB	*	GetData()	{ return	m_pData; }
	void			SetData( GenoHaploDB *pDB )	{ m_pData = pDB; }
	vector<int>		&GetNumClassN() { return m_NumClassN; }

	int		Swap( unsigned char *h1, unsigned char *h2, int nstart, int nend );
	int		Find( unsigned char *h, vector<vector<unsigned char> > &A,
			int exceptK = -1 );

	int		GetPredFreq(  unsigned char ***h );

	int		LoadData( vector<string> &filelist );
	int		LoadOutput( const char *inputfile, int nstart, int nend,
			unsigned char*** pred_h );
	int		Save( const char *inputfile, const char *outdir, int tstart, int tend );

//Added from DP
	int DeleteClass();
	int DeleteSS();
	int LoadData( const char* filename );
	void assertbug(int lineno);


	MUTDP();
	~MUTDP();


private: // for internal usage
	float	ab_h;   // = beta_h + alpha_h;
	float	a1;		// = ( alpha_h ) / (ab_h );
	float	b1;		// = ( beta_h ) / (ab_h );
	float	ab_g;	// = beta_g + alpha_g;
	float	a2;		//	= (alpha_g) / (ab_g);
	float	b2;		//	= (beta_g) / (ab_g);
	float	logB1;	// = log(B-1);
	float	logB;	// = log(B);
	float	log_mu1;	//=log(mu1);
	float	log_mu2;	//=log(mu2);	
	double	tiny;	// = pow(10, -100); //10^(-100);

	////////// hyper-parameters ///////////////
	// mutation 
	float alpha_h;
	float beta_h;

	// noisy observation
	float alpha_g;
	float beta_g;

	float mu1;
	float mu2;

	// priors for scale parameters
	float alpha_a;	//	0.5 //0.5	//1.0
	float alpha_b;	//	0.1 // 0.1	//0.4

	// only for HDP
	float gamma_a;	//	0.5	//1.0		
	float gamma_b;	//	0.1	//0.4


	// for backup
	bool	m_bClassAdded;

	int old_mk, temp_mk;
	int old_sum_mjk, temp_sum_mjk;
	int old_nk, temp_nk;
	bool old_remove_classk,temp_remove_classk;
	vector<int>	old_la0k, temp_la0k;
	vector<int>	old_la1k, temp_la1k;
	vector<int>	old_lk, temp_lk;

};

#endif // !defined(AFX_HDP_H__B8C60B76_E5B0_4C89_B1A4_86D1ABBF602C__INCLUDED_)

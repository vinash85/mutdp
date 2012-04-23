// MUTDP.cpp: implementation of the MUTDP class.
//
//////////////////////////////////////////////////////////////////////

#include <iterator>
#include <algorithm>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include "MUTDP.h"
#include "util.h"
#include "program.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////


MUTDP::MUTDP()
{
    m_pData = 0;
    m_EqClass = 0;

    //m_gamma = 0.4;

    m_doConparam = 1;		// hyper parameter sampling
    m_bCheckConvg = 1;
    m_bAmbiguous = 0;

    m_nthIter = 0;
    m_nThining = 5;
    m_nBurninIteration = 5000;
    m_nCumIteration = 5000;

    m_strOutdir = "./output/";
    m_strFilename = "test";
    m_pFPTheta = 0;

    //// parameters for prior //////
    // mutation
    alpha_h	= 80.0f;
    beta_h	= 0.8f; // 1.2f 	//	**1.2f	//* 0.8f
    // noisy observation
    alpha_g = 5000.f;
    beta_g	= 10.0f;	//	//* 5.0f
    //
    mu1		= 0.23f;
    mu2		= 0.02f;
    // scale parameter  (child DP)
    alpha_a = 0.5f; //0.5	//1.0
    alpha_b	= 0.1f; // 0.1	//0.4
    // scale parameter (master DP)
    gamma_a = 0.5f;	//1.0
    gamma_b	= 0.1f;	//0.4

}

MUTDP::~MUTDP()
{

}

int MUTDP::LoadData( vector<string> &filelist )
{
    int nfile = filelist.size();

    if ( m_pData == 0)
	m_pData = new GenoHaploDB();

    m_pData->LoadData( filelist );

    m_nBlockStart = 0;
    m_nBlockEnd = m_pData->m_numT;
    m_nBlockLength = m_nBlockEnd;

    return 1;
}


////////// initialize haplotypes and counts ///////////
int MUTDP::Init( GenoHaploDB *pDB, int nstart, int nend )
{
    int		jj;
    int		nGroups = pDB->m_nGroups;

    assert( nstart < nend );

    srand( Seed );

    // Geno-Haplo data
    m_pData = pDB;

    m_nBlockStart = nstart;
    m_nBlockEnd = nend;
    m_nBlockLength = nend - nstart;
    m_nthIter = 0;

    int		numBlockT = m_nBlockLength;

    //  constants
    ab_h = beta_h + alpha_h;
    a1 = ( alpha_h ) / (ab_h );
    b1 = ( beta_h ) / (ab_h );
    ab_g = beta_g + alpha_g;
    a2	= (alpha_g) / (ab_g);
    b2	= (beta_g) / (ab_g);
    logB1 = log(B-1);
    logB = log(B);
    log_mu1=log(mu1);
    log_mu2=log(mu2);
    tiny = pow(10, -100); //10^(-100);
    //Avinash
    clearSS();
    m_A.clear(); 

    ////////////////////////////////////////////////

    ////////////////////////////////////////////////
    // create DP: in case of known ethnic group variable
    //for ( jj = 0; jj < nGroups; jj++ )
    //{
    //DP m_dp=new DP();
    // initialize DP ethnic groups
    //for ( jj = 0; jj < nGroups; jj++ )
    //{
    //m_dp[jj].m_pDataIndex = m_pData->m_DataInGroup[jj];
    //}
    //m_pDataIndex = m_pData->m_DataInGroup;

    int I = m_pData->m_numTotalI;

    if ( m_EqClass == 0 )
    {
	m_pData->Alloc2DMemory( &m_EqClass, I, 2 );
    }

    return 1;
}


////////// initialize haplotypes and counts ///////////
int MUTDP::Initialize( GenoHaploDB *pDB, int nstart, int nend )
{
    int		ii, tt, ee, it, cc;
    int		nGroups = pDB->m_nGroups;
    bool	new_class;

    assert( nstart < nend );

    // Geno-Haplo data
    m_pData = pDB;

    m_nBlockStart = nstart;
    m_nBlockEnd = nend;
    m_nBlockLength = nend - nstart;

    int		numBlockT = m_nBlockLength;

    unsigned char ***h = m_pData->m_Haplotypes;
    unsigned char ***g = m_pData->m_Genotypes;
    int	**g_match = m_pData->m_g_match;
    int	**g_miss1 = m_pData->m_g_miss1;
    int	**g_miss2 = m_pData->m_g_miss2;
    int	***h_count = m_pData->m_h_count;
    vector<int> &n = m_NumClassN;
    //vector<int> &sum_mj = m_nSumClassN;

    //clearSS();
    //m_A.clear();

    // initialize with 1 class
    AddClass( 0, numBlockT );

    //TODO
    double tmp_gamma = m_gamma;
    m_gamma = 1;

    ////////////////////////////////////////////////
    // create DP: in case of known ethnic group variable
    //for ( jj = 0; jj < nGroups; jj++ )
    //{
    //}
    // initialize DP ethnic groups
    //for ( jj = 0; jj < nGroups; jj++ )
    //{
    //}

    //TODO
    //m_pDataIndex = m_pData->m_DataInGroup;
    int I = m_pData->m_numTotalI;

    if ( m_EqClass == 0 )
	m_pData->Alloc2DMemory( &m_EqClass, I, 2);

    // initial assignment
    for ( ii = 0; ii < I ; ii++ )
    {
	// init h
	for ( tt = nstart; tt < nend; tt++)
	{
	    float aa = rand()/(float)RAND_MAX;
	    int aa1 = 0;
	    if ( aa > 0.5 )		aa1 = 1;

	    h[0][ii][tt] = g[aa1][ii][tt];
	    h[1][ii][tt] = g[1-aa1][ii][tt];
	    if ( h[0][ii][tt] == 255 )//TODO??
	    {
		h[0][ii][tt] = 0;
		h[1][ii][tt] = 0;
	    }
	}

	// init g_match, g_miss1, g_miss2
	for (tt = nstart; tt < nend; tt++)
	{
	    int g_id0 = MIN( g[0][ii][tt], g[1][ii][tt] );
	    int g_id1 = MAX( g[0][ii][tt], g[1][ii][tt] );
	    int h_id0 = MIN( h[0][ii][tt], h[1][ii][tt] );
	    int h_id1 = MAX( h[0][ii][tt], h[1][ii][tt] );

	    if ( g_id0 == h_id0 && g_id1 == h_id1 )
	    {
		g_match[ii][tt] = 1;
		g_miss1[ii][tt] = g_miss2[ii][tt] = 0;
	    }
	    else if ( g_id0 != h_id0 && g_id1 != h_id1 )
	    {
		g_miss2[ii][tt] = 1;
		g_miss1[ii][tt] = g_match[ii][tt] = 0;
	    }
	    else
	    {
		g_miss1[ii][tt] = 1;
		g_miss2[ii][tt] = g_match[ii][tt] = 0;
	    }
	}
    }

    // calculate u[0:2] for DP
    CalNumClassU();


    // rand permutation of ii
    vector<int> iivec( I );
    for ( ii = 0; ii < I; ii++ )	iivec[ii] = ii;
    random_shuffle( iivec.begin(), iivec.end() );


    // initialization other variables
    int ** c = m_EqClass;
    for ( it = 0; it < I ; it++ )
    {
	for ( ee = 0; ee < 2; ee++ )
	{
	    ii = iivec[it];
	    int K = m_NumClassN.size();

	    //int		&m	= m_NumClassN;
	    vector<vector<int> > 	&l	= m_NumClassL;
	    vector<vector<int> >	*la = m_NumClassLA;
	    int		*u	= m_NumClassU;

	    double alpha0 = 0.7;

	    /// 1. Sample c(i,e)
	    //TODO remove FromTopLevel
	    //avinash changed Sample_EqClass_Init to Sample_EqClass
	    //cc = Sample_EqClass_Init( h[ee][ii], l, alpha0);
	    cc = Sample_EqClass( h[ee][ii], l, alpha0);

	    c[ii][ee] = cc;

	    if ( cc < K )
	    {
		//m[ cc ] ++;
		n[ cc ] ++;
		//sum_mj[ cc ] ++;
		new_class = 0;

	    }
	    else
	    {
		new_class = 1;
		// add class
		AddClass( 1, numBlockT );	// top-level
		//for ( int jt=0; jt < m_dp.size(); jt++ ) // bottom-level
		//{
		//if ( jt == jj )
		//m_dp[jt].AddClass( 1, numBlockT );
		//else
		//m_dp[jt].AddClass( 0, numBlockT );
		//}
		K++;
	    }

	    // update LA
	    for ( tt = nstart; tt < nend; tt++)
	    {
		la[ h[ee][ii][tt] ][cc][ tt-nstart ] ++;
	    }

	    m_A[cc] = Sample_A( cc, h[ee][ii], new_class );

	    K = m_A.size();

	    //
	    for ( tt = nstart; tt < nend; tt++)
	    {
		la[ h[ee][ii][tt] ][cc][ tt-nstart ] --;
	    }

	    ////////

	    // sample the haplotype "h"
	    Sample_H( h[ee][ii], h[1-ee][ii],
		    g[0][ii], g[1][ii], m_A[cc],
		    g_match[ii], g_miss1[ii], g_miss2[ii],
		    n[cc], l[cc], la[0][cc], la[1][cc],
		    h_count[ee][ii],
		    u, I );

	}	// end of e
    }	// end of i

    m_gamma = tmp_gamma;

    return 1;
}

////////////////////////////////////////////////////////////
/// do burnin iterations for Gibbs sampling combined with MH update for C variables
int MUTDP::Iterate_det_Gibbs_Met( int numIter, bool* bDone )
{
    int		ii,tt, ee, it, kk, bb;
    int		nGroups = m_pData->m_nGroups;
    bool	new_class;
    float	ErrT=0, ErrI=0, ErrSW=0;
    int		nstart = m_nBlockStart;
    int		nend = m_nBlockEnd;
    int		numBlockT = m_nBlockLength;

    // DB
    unsigned char ***h = m_pData->m_Haplotypes;
    unsigned char ***g = m_pData->m_Genotypes;
    int	**g_match = m_pData->m_g_match;
    int	**g_miss1 = m_pData->m_g_miss1;
    int	**g_miss2 = m_pData->m_g_miss2;
    int	***h_count = m_pData->m_h_count;
    vector<int> &n = m_NumClassN;
    //vector<int> &sum_mj = m_nSumClassN;

    int I = m_pData->m_numTotalI;


    m_nthIter = 0;
    printf( "##############################\n");
    printf( "Iterations for SNPs %d - %d\n", nstart, nend );
    printf( "##############################\n");
    printf("#Iter	Alpha	Gamma	K\n" );

    // iterate
    for ( int iter=0; iter < numIter; iter++ )
    {
	int K = m_NumClassN.size();

	m_Remove_Class.clear();
	for ( kk=0; kk < K; kk++)
	    m_Remove_Class.push_back(0);

	vector<bool>	&remove_class = m_Remove_Class;

	// rand permutation of ii (individuals)
	vector<int> iivec( I );
	for ( ii = 0; ii < I; ii++ )	iivec[ii] = ii;
	random_shuffle( iivec.begin(), iivec.end() );

	// initialization other variables
	int ** c = m_EqClass;
	for ( it = 0; it < I ; it++ )
	{
	    for ( ee = 0; ee < 2; ee++ )
	    {
		ii = iivec[it];

		if ( !(bDone && bDone[ii]) )
		{
		    K = m_NumClassN.size();

		    // Ethinic Group selection
		    //jj = m_pData->m_EthnicGroup[ii];

		    //vector<int>				&m	= m_dp[jj].m_NumClassN;
		    vector<vector<int> >	&l	= m_NumClassL;
		    vector<vector<int> >	*la = m_NumClassLA;
		    int						*u	= m_NumClassU;

		    double alpha0 = m_alpha; //TODO is 1 in DP constructor

		    int	old_c = c[ii][ee];
		    // in BackOldSS variables n, la, l are backed up


		    BackupOldSS( old_c, remove_class );
		    // in DeleteSS variables n[cc]--, la[h]--, l are backed up. This is done for gibbs sampling because for sampling of cc the contribution due to cc needs to be removed.

		    DeleteSS( ii, ee, c[ii][ee]);

		    old_lk = l[ old_c ];

		    if ( n[ old_c ] == 0 )
			remove_class[ old_c ] = 1;
		    else
			remove_class[ old_c ] = 0;

		    //////////////////////////////////
		    /// 1. Sample c(i,e)

		    int cc = Sample_EqClass( h[ee][ii], l, alpha0);

		    c[ii][ee] = cc;

		    // for back-up
		    if ( c[ii][ee] < K )
		    {
			BackupUpdateSS( cc, remove_class, 0 );

			new_class = 0;
			remove_class[ cc ] = 0;

			n[ cc ] ++;
			//sum_mj[ cc ] ++;



			for ( tt = nstart; tt < nend; tt++)
			{
			    la[ h[ee][ii][tt] ][cc][tt-nstart]++;
			}
		    }
		    else
		    {
			new_class = 1;
			kk = 0;
			while ( remove_class[kk] == 0 && kk < K )
			{
			    kk++;
			}
			if ( kk < K )
			{
			    cc = kk;

			    BackupUpdateSS( cc, remove_class, 0 );

			    c[ii][ee] = cc;
			    remove_class[ cc  ] = 0;
			    n[ cc ] = 1;
			    //sum_mj[ cc ] = 1;
			    n[ cc ] = 1;

			    for ( tt=0; tt < numBlockT; tt++)
			    {
				int hh = h[ee][ii][tt];
				l[ cc ][tt] = 0;
				la[ hh ][ cc ][tt] = 1;
				la[1-hh][ cc ][tt] = 0;
			    }
			}
			else
			{
			    BackupUpdateSS( cc,  remove_class, 1 );

			    // add class
			    AddClass( 1, numBlockT );	// top-level


			    remove_class.push_back(0);

			    for ( tt = nstart; tt < nend; tt++ )
			    {
				la[ h[ee][ii][tt] ][cc][ tt-nstart ] ++;
			    }

			    K++;
			}
		    }

		    vector<unsigned char> old_a = m_A[old_c];

		    // sample ancestor A
		    vector<unsigned char> temp_a = Sample_A( cc, h[ee][ii], new_class );

		    // Metropolis test 
		    bool test = TestAcceptance( old_c, cc, h[ee][ii], old_a, temp_a, n, l );

		    if ( test )
		    {
			m_A[cc] = temp_a;
		    }
		    else
		    {
			c[ii][ee] = old_c;
			RollBack( ii, ee, old_c, cc );
		    }

		    //
		    cc = c[ii][ee];
		    for (bb = 0; bb < B; bb++)
		    {
			for ( tt = nstart; tt < nend; tt++)
			{
			    if (  bb == h[ee][ii][tt] )
				la[bb][cc][tt - nstart]--;
			}
		    }

		    // sample the haplotype "h"
		    Sample_H( h[ee][ii], h[1-ee][ii],
			    g[0][ii], g[1][ii], m_A[cc],
			    g_match[ii], g_miss1[ii], g_miss2[ii],
			    n[cc], l[cc], la[0][cc], la[1][cc],
			    h_count[ee][ii],
			    u, I);
		}
	    }	// end of e
	}	// end of i
	// Date 20 april
	CondenseList( 1 );

	K = m_NumClassN.size();

	//	if ( iter % 5 == 0 )
	{
	    MergeList();
	    CondenseList();//TODO
	}

	if ( m_pFPTheta)
	    EstimateTheta( iter );

	if ( m_doConparam && iter >= 100 )
	    Sample_Conparam( 1, 10 );

	if ( iter % 100 == 0 )
	{
	    printf("%d	%.4f %.4f	%d\n", m_nthIter, m_alpha, m_gamma,  K );
	}

	m_nthIter++;

    }	// end of iter

    return 1;
}


int MUTDP::Iterate_cum_Gibbs_Met( int numIter, bool* bDone )
{
    int		ii, tt, ee, it, kk;
    int		nGroups = m_pData->m_nGroups;
    bool	new_class;

    int		nstart = m_nBlockStart;
    int		nend = m_nBlockEnd;
    int		numBlockT = m_nBlockLength;

    // DB
    unsigned char ***h = m_pData->m_Haplotypes;
    unsigned char ***g = m_pData->m_Genotypes;
    int	**g_match = m_pData->m_g_match;
    int	**g_miss1 = m_pData->m_g_miss1;
    int	**g_miss2 = m_pData->m_g_miss2;
    int	***h_count = m_pData->m_h_count;
    vector<int> &n = m_NumClassN;
    //vector<int> &sum_mj = m_nSumClassN;

    int I = m_pData->m_numTotalI;

    // initialize cum_h
    int *cross = new int[I];
    memset( cross, 0, sizeof(int)*I );
    m_pData->UpdateCumHaplotype( h, 0, nstart, nend, &cross );

    // backup h_old
    unsigned char **h_old[2];
    m_pData->Alloc2DMemory( &(h_old[0]), I, numBlockT );
    m_pData->Alloc2DMemory( &(h_old[1]), I, numBlockT );

    for (ii=0; ii<I; ii++)
    {
	memcpy( h_old[0][ii], &(h[0][ii][nstart]),
		sizeof(unsigned char)*numBlockT );
	memcpy( h_old[1][ii], &(h[1][ii][nstart]),
		sizeof(unsigned char)*numBlockT );
    }

    unsigned char **predh_old[2];
    if ( m_bCheckConvg )
    {
	m_pData->Alloc2DMemory( &(predh_old[0]), I, numBlockT );
	m_pData->Alloc2DMemory( &(predh_old[1]), I, numBlockT );

	for (ii=0; ii<I; ii++)
	{
	    memcpy( predh_old[0][ii], &(h[0][ii][nstart]),
		    sizeof(unsigned char)*numBlockT );
	    memcpy( predh_old[1][ii], &(h[1][ii][nstart]),
		    sizeof(unsigned char)*numBlockT );
	}

    }
    vector<int> &traceDiff = m_traceDiff;

    // iterate
    for ( int iter=0; iter < numIter; iter++ )
    {
	int K = m_NumClassN.size();

	m_Remove_Class.clear();
	for ( kk=0; kk < K; kk++)
	    m_Remove_Class.push_back(0);

	vector<bool>	&remove_class = m_Remove_Class;

	///////////////////////////
	// Random permutation of individuals
	vector<int> iivec( I );
	for ( ii = 0; ii < I; ii++ )	iivec[ii] = ii;
	random_shuffle( iivec.begin(), iivec.end() );

	//////////////////////////////////
	// Main iteration
	int ** c = m_EqClass;
	for ( it = 0; it < I ; it++ )
	{
	    for ( ee = 0; ee < 2; ee++ )
	    {
		ii = iivec[it];

		if ( !(bDone && bDone[ii]) )
		{
		    K = m_NumClassN.size();

		    // Ethinic Group selection
		    //jj = m_pData->m_EthnicGroup[ii];

		    vector<int>				&m	= m_NumClassN;
		    vector<vector<int> >	&l	= m_NumClassL;
		    vector<vector<int> >	*la = m_NumClassLA;
		    int	*u = m_NumClassU;

		    double alpha0 = m_alpha;

		    int	old_c = c[ii][ee];

		    //////////////////////////////////
		    // for roll-back in MH updating
		    BackupOldSS( old_c, remove_class );
		    // for gibbs sampling n_k and m_jk are updated for sampling of c_ie; its contribution is removed
		    DeleteSS( ii, ee, c[ii][ee] );
		    old_lk = l[ old_c ];

		    if ( m[ old_c ] == 0 )
			remove_class[ old_c ] = 1;
		    else
			remove_class[ old_c ] = 0;

		    //////////////////////////////////
		    /// 1. Sample c(i,e)
		    //cc: current class
		    //ee: mother fater
		    //ii : individual
		    //tt: block start to block end
		    //kk : ancestor
		    //c(i,e): equivalence class. h(i) =a(c(i,e)) + noise
		    //jj: ethnic group
		    int cc = Sample_EqClass( h[ee][ii], l, alpha0);
		    c[ii][ee] = cc;

		    // update

		    if ( c[ii][ee] < K ) //old color
		    {
			BackupUpdateSS( cc, remove_class, 0 );
			new_class = 0;
			remove_class[ cc ] = 0;

			n[ cc ] ++;

			//n[cc] = n_k

			for ( tt = nstart; tt < nend; tt++)
			    la[ h[ee][ii][tt] ][cc][tt-nstart]++;
		    }
		    else
		    {
			new_class = 1;
			kk = 0;
			while ( remove_class[kk] == 0 && kk < K )
			{
			    kk++;
			}
			// if there is a removed class, put the new class into it.
			if ( kk < K )
			{
			    cc = kk;

			    BackupUpdateSS( cc, remove_class, 0 );

			    c[ii][ee] = cc;
			    remove_class[ cc  ] = 0;
			    //m[ cc ] = 1;
			    //sum_mj[ cc ] = 1;
			    n[ cc ] = 1;

			    for ( tt=0; tt < numBlockT; tt++)
			    {
				int hh = h[ee][ii][tt];
				l[ cc ][tt] = 0; // l(k,t) = 0 -->  a(k,t) = h(k,t)
				la[ hh ][ cc ][tt] = 1;// la
				la[1-hh][ cc ][tt] = 0;
			    }
			}
			else
			{
			    BackupUpdateSS( cc, remove_class, 1 );

			    // add class
			    AddClass( 1, numBlockT );	// top-level

			    //for ( int jt=0; jt < m_dp.size(); jt++ ) // bottom-level
			    //{
			    //if ( jt == jj )
			    //m_dp[jt].AddClass( 1, numBlockT );
			    //else
			    //m_dp[jt].AddClass( 0, numBlockT );
			    //}

			    remove_class.push_back(0);

			    // update SS
			    for ( tt = nstart; tt < nend; tt++ )
				la[ h[ee][ii][tt] ][cc][ tt-nstart ] ++;

			    K++;
			}
		    }

		    vector<unsigned char> old_a = m_A[old_c];

		    /////////////////////////////
		    // 2. Sample ancestor Ak
		    vector<unsigned char> temp_a = Sample_A( cc, h[ee][ii], new_class );

		    // Metropolis-Hasting test
		    if ( TestAcceptance( old_c, cc, h[ee][ii], old_a, temp_a, n, l ) )
		    {
			m_A[cc] = temp_a;
		    }
		    else
		    {
			c[ii][ee] = old_c;
			RollBack( ii, ee, old_c, cc );
		    }

		    //
		    cc = c[ii][ee];
		    for ( tt = nstart; tt < nend; tt++)
		    {
			la[ h[ee][ii][tt] ][cc][tt - nstart]--;
		    }

		    ///////////////////////////
		    // 3. sample Haplotype H(i,e)
		    Sample_H( h[ee][ii], h[1-ee][ii],  g[0][ii], g[1][ii], m_A[cc],
			    g_match[ii], g_miss1[ii], g_miss2[ii],
			    m[cc], l[cc], la[0][cc], la[1][cc],
			    h_count[ee][ii],
			    u, I);
		}
	    }	// end of e
	}	// end of i


	CondenseList( 1 );		// eliminate duplication

	K = m_NumClassN.size();

	if ( iter % 100 == 0 )
	{
	    printf("%d	%.4f %.4f	%d\n", m_nthIter, m_alpha, m_gamma,  K );
	}

	//	if ( iter % 5 == 0 )
	{
	    MergeList();
	    CondenseList();
	}

	if ( m_pFPTheta)
	    EstimateTheta( iter );

	if ( m_doConparam )
	    Sample_Conparam( 1, 10 );

	// add samples for posterior mean
	if ( iter % m_nThining == 0 )
	{
	    m_pData->UpdateCumHaplotype( h, h_old, nstart, nend, &cross );
	}

	//
	for (ii=0; ii<I; ii++)
	{
	    memcpy( h_old[0][ii], &(h[0][ii][nstart]),
		    sizeof(unsigned char)*numBlockT );
	    memcpy( h_old[1][ii], &(h[1][ii][nstart]),
		    sizeof(unsigned char)*numBlockT );
	}

	m_nthIter++;
    }

    // delete memory
    for (ii=0 ; ii < I; ii++)
    {
	if ( h_old[0][ii] ) delete [] h_old[0][ii];
	if ( h_old[1][ii] ) delete [] h_old[1][ii];
    }
    if ( h_old[0] ) delete [] h_old[0];
    if ( h_old[1] ) delete [] h_old[1];

    if ( cross ) delete [] cross;

    return 1;
}

//////////////////////////////////////////////////////////
// If there are duplications in Ancestoral pool, 
// merge them into one and update the corresponding counters
///////////////////////////////////////////////////////////
int MUTDP::MergeList()
{
    int K = m_A.size();
    int kk, tt, ka;

    vector<vector<unsigned char> > A_new;
    vector<int> MapK(K,0);
    int ind = 0;
    vector<bool> bMapped(K,0);

    for (kk=0; kk<K; kk++)
    {
	if ( bMapped[kk] == 0 )
	{
	    m_Remove_Class[kk] = 0;
	    MapK[kk] = kk;
	    bMapped[kk] = 1;
	    vector<unsigned char> cmp( m_A[kk] );
	    for (ka=kk+1; ka<K; ka++)
	    {
		if ( bMapped[ka] == 0 )
		{
		    if ( m_A[ka] == m_A[kk] )
		    {
			m_Remove_Class[ka] = 1;
			MapK[ka] = kk;
			bMapped[ka] = 1;

			m_NumClassN[kk] += m_NumClassN[ka];
			//m_nSumClassN[kk] += m_nSumClassN[ka];
			for (tt=0; tt<m_nBlockLength; tt++)
			{
			    m_NumClassLA[0][kk][tt] += m_NumClassLA[0][ka][tt];
			    m_NumClassLA[1][kk][tt] += m_NumClassLA[1][ka][tt];
			    m_NumClassL[kk][tt] += m_NumClassL[ka][tt];
			}
		    }
		}
	    }
	    ind++;
	}
    }

    for (int ii=0; ii < m_pData->m_numTotalI; ii++)
    {
	for( int ee=0; ee<2; ee++)
	{
	    int cc = m_EqClass[ii][ee] ;
	    assert( MapK[cc] >= 0 );
	    m_EqClass[ii][ee] = MapK[cc];

	}

	return 1;
    }
}


///////////////////////////////////////////
// Remove the empty class from the list
//	and update the counters and variables
int MUTDP::CondenseList( int bBackupA )
{
    int ii, ee, kk;
    int K = m_NumClassN.size();

    int ind = 0, rem_k;
    vector<int>	MapK(K,0);
    bool foundR = 0;
    kk = 0;
    while ( kk < K )
    {
	if ( m_Remove_Class[kk] == 0 )
	    //		if ( m_nSumClassN[kk] > 0 )
	{
	    MapK[kk] = ind;
	    ind++;
	}
	else
	{
	    MapK[kk] = -1;
	    if ( foundR == 0 )
	    {
		foundR = 1;
		rem_k = ind;
	    }
	}
	kk++;
    }

    if ( foundR == 0 )
	return 1;


    ///////////
    for (kk=K-1; kk>=rem_k; kk--)
    {
	if ( MapK[kk] < 0 )
	{
	    if (  bBackupA == 1 )
		m_Buffer_A.push_back( m_A[kk] );

	    m_NumClassN.erase( m_NumClassN.begin()+kk );
	    //m_nSumClassN.erase( m_nSumClassN.begin()+kk );
	    m_A.erase( m_A.begin()+kk );

	    //			m_Remove_Class.erase( m_Remove_Class.begin()+kk );

	    //for (jj=0; jj<J; jj++)
	    //{
	    m_NumClassN.erase( m_NumClassN.begin()+kk );
	    m_NumClassL.erase( m_NumClassL.begin()+kk );
	    m_NumClassLA[0].erase( m_NumClassLA[0].begin()+kk );
	    m_NumClassLA[1].erase( m_NumClassLA[1].begin()+kk );
	    //}

	}
    }

    int I = m_pData->m_numTotalI;

    for ( ii=0; ii < I; ii++)
    {
	for ( ee=0; ee<2; ee++)
	{
	    int cc = m_EqClass[ii][ee];
	    assert( MapK[cc] >= 0 );
	    m_EqClass[ii][ee] = MapK[cc];
	}
    }

    return 1;
}

//////////////////////////////////////////////////////
//	Sample H
//		h: haplotype to be sampled
//		h1: complementary haplotype to h
//		g0, g1: genotypes
//		ak: ancestral haplotype for h to be inherited
//		g_match, g_miss1, g_miss2: ss for H->G observation
//		mk, lk, lak0, lak1, h_count, u: ss
//		Ij: the number of individuals of group j
///////////////////////////////////////////////////////
int MUTDP::Sample_H( unsigned char *h, unsigned char *h1,
	unsigned char *g0, unsigned char *g1,
	vector<unsigned char>  &ak,
	int	*g_match, int *g_miss1, int *g_miss2,
	int nk, vector<int> &lk, vector<int> &lak0, vector<int> &lak1,
	int *h_count,
	int *u, int I)
{
    int tt, bb;
    int numBlockT = m_nBlockLength;
    int nstart = m_nBlockStart;
    int	nend = m_nBlockEnd;
    int minh, maxh, ming, maxg;

    vector<int> &n = m_NumClassN;
    int K = n.size();

    double log_a2 = log(a2);
    double log_b2 = log(b2);
    vector<double> log_pH;

    for ( tt = nstart; tt < nend; tt++ )
    {
	int u_temp[3][B];
	memset( (int*)u_temp, 0, sizeof(int)*3*B );
	u[0] = u[0] - g_match[tt];
	u[1] = u[1] - g_miss1[tt];
	u[2] = u[2] - g_miss2[tt];
	for (bb = 0; bb < B; bb++)
	{
	    minh = MIN( bb, h1[tt] );
	    maxh = MAX( bb, h1[tt] );
	    ming = MIN( g0[tt], g1[tt] );
	    maxg = MAX( g0[tt], g1[tt] );
	    if ( (minh == ming) && (maxh == maxg) )
		u_temp[0][bb] = 1;
	    else if ( ( minh != ming ) && (maxh != maxg ) )
		u_temp[2][bb] = 1;
	    else
		u_temp[1][bb] = 1;
	}
	int l_temp[] = {0, 0};
	l_temp[ ak[tt-nstart] ] = 1;

	log_pH.clear();
	double min_log_pH = INF;
	for (bb = 0; bb < B; bb++)
	{
	    double lp = lgamma( alpha_h + lk[tt-nstart] + l_temp[bb] )
		+ lgamma( beta_h + nk - lk[tt-nstart] - l_temp[bb])
		- lgamma( nk + ab_h )
		- ( nk - lk[tt-nstart] - l_temp[bb])*logB1 +
		( u[0] + u_temp[0][bb] )*log_a2 +
		( I*m_nBlockLength - u[0] - u_temp[0][bb])*log_b2
		+ ( u[1] + u_temp[1][bb] )*log_mu1
		+ ( u[2] + u_temp[2][bb] )*log_mu2;

	    log_pH.push_back( lp );
	    if ( lp < min_log_pH )
		min_log_pH = lp;
	}
	for (bb=0; bb<B; bb++)
	    log_pH[bb] = exp( log_pH[bb] - min_log_pH );

	normalize( log_pH );
	// sample ht
	h[tt] = sample_discrete( log_pH );

	if ( ak[tt-nstart] == h[tt] )
	    h_count[tt] = 1;
	else
	    h_count[tt] = 0;

	bb = h[tt];
	minh = MIN( bb, h1[tt] );
	maxh = MAX( bb, h1[tt] );
	ming = MIN( g0[tt], g1[tt] );
	maxg = MAX( g0[tt], g1[tt] );
	if ( (minh == ming) && (maxh == maxg) )
	{
	    g_match[tt] = 1;
	    g_miss1[tt] = 0;
	    g_miss2[tt] = 0;

	    u[0] += 1;
	}
	else if ( ( minh != ming ) && (maxh != maxg ) )
	{
	    g_match[tt] = 0;
	    g_miss1[tt] = 0;
	    g_miss2[tt] = 1;

	    u[2] += 1;
	}
	else
	{
	    g_match[tt] = 0;
	    g_miss1[tt] = 1;
	    g_miss2[tt] = 0;

	    u[1] += 1;
	}

	lk[tt-nstart] += h_count[tt];

	if ( h[tt] == 0 )
	    lak0[tt-nstart] += 1;
	else
	    lak1[tt-nstart] += 1;

    }

    return 1;
}


//new_class is flag
vector<unsigned char> MUTDP::Sample_A( int cc, unsigned char *h, bool new_class )
{
    int tt, bb;
    int numBlockT = m_nBlockLength;
    vector<unsigned char> ak;
    int at;

    // previously instantiated ancestor
    if ( !new_class )
    {
	for ( tt = 0; tt < numBlockT; tt++ )
	{
	    double LA_all[2], log_pA;
	    vector<double> pA;
	    for ( bb = 0; bb < B; bb++ )
	    {
		//
		LA_all[bb] = m_NumClassLA[bb][cc][tt];

		double LA_1 = m_NumClassN[cc] - LA_all[bb];
		//eqn:8

		log_pA = lgamma( alpha_h + LA_all[bb] ) +
		    lgamma( beta_h + LA_1 )
		    - lgamma ( m_NumClassN[cc] + ab_h ) - LA_1*logB1;

		pA.push_back( exp( log_pA ) );
	    }
	    normalize( pA );
	    at = sample_discrete( pA );
	    if ( at < 0 && at > B )
		assert(0);
	    ak.push_back( (unsigned char)at );
	}
    }
    else		// newly added class
    {
	vector<double> pA, ppA;
	double pt = beta_h / ( (B-1)*(ab_h) );
	for ( bb = 0; bb < B; bb++)
	    ppA.push_back( pt );

	for ( tt = 0; tt < numBlockT; tt++ )
	{
	    pA.clear();
	    pA = ppA;
	    if ( h[tt+m_nBlockStart] != B )
		pA[ h[tt+m_nBlockStart] ] = alpha_h / ab_h;

	    normalize( pA );
	    at = sample_discrete( pA ) ;
	    if ( at < 0 && at > B )
		assert(0);
	    ak.push_back( (unsigned char)at );
	}

    }

    return ak;
}


int MUTDP::CalNumClassU()
{
    int ii, tt;
    unsigned char ***h = m_pData->m_Haplotypes;
    unsigned char ***g = m_pData->m_Genotypes;
    int		**g_match = m_pData->m_g_match;
    int		**g_miss1 = m_pData->m_g_miss1;
    int		**g_miss2 = m_pData->m_g_miss2;

    int nGroups = m_pData->m_nGroups;

    int		*uj		= m_NumClassU;
    int		I	= m_pData->m_numTotalI;
    //int	&pIndex = m_pDataIndex;
    uj[0] = uj[1] = uj[2] = 0;

    for ( int ii = 0; ii < I; ii++ )
    {
	    //ii = pIndex[it];
	for ( tt = m_nBlockStart; tt < m_nBlockEnd; tt++ )
	{
	    if ( g_match[ii][tt] == 1 )
		uj[0] ++;
	    else if ( g_miss1[ii][tt] == 1 )
		uj[1] ++;
	    else if ( g_miss2[ii][tt] == 1 )
		uj[2] ++;
	}
    }
}

///////////////////////////////////////////////////////
//		Propose C(i,e) from prior(tau)
//		m = DP[j].numClassN (number of haplotypes of  Ancestor in each ethnic group)
//		l = numClassL
//		FromTopLevel = whether the draw is from top or bottom urn
//		@returns: which color is drawn Cie
///////////////////////////////////////////////////////
//int MUTDP::Sample_EqClass( unsigned char *h, vector<int> &m,
int MUTDP::Sample_EqClass_Init( unsigned char *h, 
	vector<vector<int> > &l, double alpha0)
{
    int cc, kk;
    // n = number of haplotypes for a given ancestor
    vector<int> &n = m_NumClassN;
    int K = n.size();
    //double a_nsumg = alpha0 / ( sum( n ) + m_gamma ) ;

    // DP_vec stores discrete prob dist of phi's (each \phi denotes random variable associated with the color varying from 1..K+1, K+1 is new color, it is associated with lower urn) //sum(DP_Vec) = 1
    // alpha0 = tau
    vector<double>	DP_vec;

    //this is proababilty of choosing old color.
    for (kk=0; kk < K; kk++)
    {
	DP_vec.push_back( n[kk] );
    }
    //this is proababilty of choosing new color //eqn:2
    DP_vec.push_back( alpha0 );

    normalize( DP_vec );

    cc = sample_discrete( DP_vec );

    // determine whether from bottom level or from top level

    return cc;
}


//int MUTDP::Sample_EqClass_Init(unsigned char *h,
int MUTDP::Sample_EqClass(unsigned char *h,
	vector<vector<int> > &l, double alpha0)
{
    int cc;
    int kk, tt;
    vector<int> &n = m_NumClassN;
    int K = n.size();


    // ph | a
    double *log_ph = new double[K+1];
    double *log_ph_temp = new double[K+1];

    for ( kk = 0 ; kk < K ; kk++ )
    {
	double	mc = log( n[kk] + ab_h );
	double	mc1 = mc + logB1;
	log_ph[kk] = 0;
	if ( n[kk] != 0 )
	{
	    double log_pht;
	    for (tt = m_nBlockStart; tt < m_nBlockEnd; tt++)
	    {
		if ( m_A[kk][tt - m_nBlockStart] == h[tt] )
		    log_pht = log( alpha_h + l[kk][tt-m_nBlockStart] ) - mc;
		else
		    log_pht = log( beta_h + n[kk] - l[kk][tt-m_nBlockStart] ) - mc1;
		log_ph[kk] += log_pht;
	    }
	    log_ph_temp[kk] = log_ph[kk];
	}
	else
	{
	    log_ph[kk] = -INF;
	    log_ph_temp[kk] = INF;
	}
    }
    log_ph[K] = -m_nBlockLength * logB;
    log_ph_temp[K] = log_ph[K];

    double *log_DP_predict = new double[K+1];
    double *log_DP_predict_temp = new double[K+1];
    double min_DP_predict = INF;

    // DP_vec
    vector<double>	DP_vec;
    DP_vec.assign( K+1, 0 );
    int sumn = (int)sum<int>( n );

    for ( kk = 0; kk < K; kk++)
    {
	//DP_vec[kk] = m[kk] + alpha0 * n[kk]/(sumn + m_gamma);
	DP_vec[kk] = n[kk];
    }
    DP_vec[K] = alpha0;

    normalize( DP_vec );

    for ( kk = 0; kk < K+1; kk++ )
    {
	double logdp		= log( DP_vec[kk] + tiny );

	log_DP_predict[kk]	= log_ph[kk] + logdp;
	log_DP_predict_temp[kk] = log_ph_temp[kk] + logdp;

	if ( log_DP_predict_temp[kk] < min_DP_predict )
	    min_DP_predict = log_DP_predict_temp[kk];
    }

    vector<double> DP_predict;
    for ( kk = 0; kk < K+1; kk++ )
    {
	double tmp = log_DP_predict[kk] - min_DP_predict;
	if ( tmp < -30 ) tmp = 0;
	else	tmp = exp( tmp );

	DP_predict.push_back( tmp  );
    }
    normalize( DP_predict );

    cc = sample_discrete( DP_predict );


    if ( log_ph ) delete [] log_ph;
    if ( log_ph_temp ) delete [] log_ph_temp;
    if ( log_DP_predict ) delete [] log_DP_predict;
    if ( log_DP_predict_temp ) delete [] log_DP_predict_temp;

    return cc;
}

int MUTDP::clearSS()
{
    m_NumClassN.clear();

    return 1;
}

int MUTDP::AddClass( int initvalue, int numT )
{
    // initialize with 1 class
    m_NumClassN.push_back( initvalue );

    vector<unsigned char> zeros1( numT, 0 );
    m_A.push_back( zeros1 );
    vector<int> zeros( numT, 0 );
    m_NumClassL.push_back( zeros );
    m_NumClassLA[0].push_back( zeros );
    m_NumClassLA[1].push_back( zeros );

    return 1;
}

int MUTDP::DeleteSS( int ii, int ee, int cc)
{
    unsigned char *h = m_pData->m_Haplotypes[ee][ii];
    int *h_count = m_pData->m_h_count[ee][ii];

    //global m_NumClassN is n_k , number of colors of ancestor

    //m_nSumClassN[ cc ] -= 1;
    m_NumClassN[ cc ] -= 1;


    for (int tt = 0; tt < m_nBlockLength; tt++)
    {
	m_NumClassLA[ h[tt+m_nBlockStart]][cc][tt] -= 1;//m_NumClassLA = m_jk//TODO
	m_NumClassL[ cc ][tt] -= h_count[ tt+m_nBlockStart ];
    }

    return 1;
}

bool MUTDP::TestAcceptance( int old_c, int new_c, unsigned char *h,
	vector<unsigned char> &old_a, vector<unsigned char> &temp_a,
	vector<int> &n, vector<vector<int> > &l )
{
    if ( (old_c == new_c) && (temp_a == old_a) )
	return 1;

    //if ( old_c != new_c )
    //int test = 1;

    double	mc = log( n[old_c] + ab_h );
    double	mc1 = mc + logB1;

    double	log_ph_old = 0;
    double	log_ph_new = 0;
// m = number of descendent from each ancestor K
    double	mc_n = log( n[new_c] + ab_h );
    double	mc1_n = mc_n + logB1;

    // posterior likelihood of h(:,i,e)
    for (int tt = 0; tt < m_nBlockLength; tt++)
    {
	assert( n[old_c] >= l[old_c][tt] );

	if ( old_a[tt] == h[tt+m_nBlockStart] )
	    log_ph_old += log( alpha_h + l[old_c][tt] ) - mc;
	else
	    log_ph_old += log( beta_h + n[ old_c ] - l[old_c][tt] ) - mc1;

	assert( n[new_c] >= l[new_c][tt] );
	if ( temp_a[tt] == h[tt+m_nBlockStart] )
	    log_ph_new += log( alpha_h + l[old_c][tt] ) - mc;
	else
	    log_ph_new += log( beta_h + n[ old_c ] - l[old_c][tt] ) - mc1;
    }


    double acceptance_p = MIN( 1, exp(log_ph_new - log_ph_old) );

    if ( (double)rand()/(double)(RAND_MAX) <= acceptance_p )
	return 1;
    else
	return 0;
}

int	MUTDP::BackupOldSS( int old_c, vector<bool> &remove_class  )
{

    old_nk = m_NumClassN[ old_c ];
    //old_sum_mjk = m_nSumClassN[ old_c ];

    old_remove_classk = remove_class[ old_c ];

    //old_mk = m_NumClassN[ old_c ];

    old_la0k = m_NumClassLA[0][ old_c ];
    old_la1k = m_NumClassLA[1][ old_c ];

    old_lk = m_NumClassL[ old_c ];

    return 1;
}

// new_c: current color sampled
// classAdded : whether new color has been sampled or not
//
int	MUTDP::BackupUpdateSS( int new_c, vector<bool> &remove_class, bool classAdded  )
{
    m_bClassAdded = classAdded;

    if ( !m_bClassAdded )
    {

	temp_nk = m_NumClassN[ new_c ];

	temp_remove_classk = remove_class[ new_c ];


	temp_la0k = m_NumClassLA[0][ new_c ];
	temp_la1k = m_NumClassLA[1][ new_c ];

	temp_lk = m_NumClassL[ new_c ];
    }

    return 1;
}

int MUTDP::RollBack( int ii, int ee, int old_c, int new_c )
{
    m_NumClassN[old_c] = old_nk;
    m_Remove_Class[old_c] = old_remove_classk;

    m_NumClassLA[0][ old_c ] = old_la0k;
    m_NumClassLA[1][ old_c ] = old_la1k;
    m_NumClassL[ old_c ] = old_lk;

    if ( m_bClassAdded )
    {
	m_NumClassN.pop_back();
	m_Remove_Class.pop_back();
	m_A.pop_back();

    }
    else if ( old_c != new_c )
    {
	m_NumClassN[new_c] = temp_nk;
	m_Remove_Class[new_c] = temp_remove_classk;

	m_NumClassLA[0][ new_c ] = temp_la0k;
	m_NumClassLA[1][ new_c ] = temp_la1k;
	m_NumClassL[ new_c ] = temp_lk;
    }



    return 1;
}

////////////////////////////////////////////////////////
//	gather posterior samples and save posterior mean
////////////////////////////////////////////////////////
int MUTDP::Sample_Pred()
{
    int I = m_pData->m_numTotalI;

    int nstart = m_nBlockStart;
    int nend = m_nBlockEnd;
    //int J = m_dp.size();

    int ii, tt;

    unsigned char **pred_h[2];

    pred_h[0] = m_pData->m_Pred_Haplotypes[0];
    pred_h[1] = m_pData->m_Pred_Haplotypes[1];

    unsigned char ***h = m_pData->m_Haplotypes;
    unsigned char ***h0 = m_pData->m_TrueHaplotypes;
    unsigned char ***g = m_pData->m_Genotypes;
    unsigned char **g_raw = m_pData->m_RawGenotypes;
    int		 **cum_h[2][2];

    for( int bb=0; bb < B; bb++)
    {
	cum_h[bb][0] = m_pData->m_Cum_Haplotypes[bb][0];
	cum_h[bb][1] = m_pData->m_Cum_Haplotypes[bb][1];
    }

    for ( ii = 0 ; ii < I; ii ++ )
    {
	for ( tt=nstart; tt<nend; tt++)
	{
	    if ( cum_h[0][0][ii][tt] > cum_h[1][0][ii][tt] )
		pred_h[0][ii][tt] = 0;
	    else
		pred_h[0][ii][tt] = 1;

	    if ( cum_h[0][1][ii][tt] > cum_h[1][1][ii][tt] )
		pred_h[1][ii][tt] = 0;
	    else
		pred_h[1][ii][tt] = 1;
	}
    }

    int ct3=0, cind = 0, mis1=0;
    int hetI = 0;
    for (ii=0 ; ii < I; ii++)
    {
	int ct1 = 0;
	int ct2 = 0;
	int hetero = 0;

	for (tt = nstart ; tt < nend ; tt++)
	{
	    if ( g_raw[ii][tt] == 2 )
	    {
		unsigned char h0 = pred_h[0][ii][tt],
			      h1 = pred_h[1][ii][tt];
		if ( pred_h[0][ii][tt] == g[0][ii][tt] )
		    h1 = g[1][ii][tt];

		else if ( pred_h[0][ii][tt] == g[1][ii][tt] )
		    h1 = g[0][ii][tt];

		if ( pred_h[1][ii][tt] == g[0][ii][tt] )
		    h0 = g[1][ii][tt];

		else if ( pred_h[1][ii][tt] == g[1][ii][tt] )
		    h0 = g[0][ii][tt];

		pred_h[0][ii][tt] = h0;
		pred_h[1][ii][tt] = h1;
	    }
	    if ( g_raw[ii][tt] == 0 || g_raw[ii][tt] == 1 )
	    {
		pred_h[0][ii][tt] = g[0][ii][tt];
		pred_h[1][ii][tt] = g[1][ii][tt];
	    }
	}
    }

    return 1;
}

int MUTDP::LoadOutput( const char *inputfile, int nstart, int nend,
	unsigned char*** pred_h )
{
    int tt, ee;

    FILE *fp = fopen( inputfile, "r" );

    // 0. Load Header
    char tmp[500];
    int pos = -1;
    while ( !feof( fp ) )
    {
	fgets( tmp, 500, fp );
	string oline( tmp );
	pos = oline.find( "[H]" );
	if ( pos >= 0 ) break;
    }

    if ( feof( fp ) )
	return 0;

    int i = -1;
    pos = 0;
    while (  pos >= 0 && ~feof( fp ) )
    {
	fgets( tmp, 500, fp );
	string oline( tmp );
	pos = oline.find(  "---" );
	int post = oline.find(  "," );
	if ( pos < 0 )
	    break;

	i++;
	string stmp = oline.substr(0, post);
	int it;
	sscanf( stmp.c_str(), "%d", &it );
	//		printf( "%d, %d, %s \n", i,it, tmp );
	if ( it != i )
	    break;
	char h_ie[5000];

	for ( ee = 0; ee < 2; ee++ )
	{
	    fgets( h_ie, 5000, fp );

	    for ( tt = nstart; tt < nend; tt++ )
	    {
		char ochar = h_ie[tt-nstart];
		//				pred_h[ee][i][tt] = atoi( &(ochar) );
		if ( ochar == '0' )
		    pred_h[ee][i][tt] = 0;
		else
		    pred_h[ee][i][tt] = 1;
		////////////////////////////////
		m_pData->m_Pred_Haplotypes[ee][i][tt] = pred_h[ee][i][tt] ;
		////////////////////////////////
	    }
	}
	if ( feof(fp) )
	    break;
    }


    GetPredFreq( m_pData->m_Pred_Haplotypes  );
    fclose(fp);

    return 1;
}

/*
   int MUTDP::Save( const char *inputfile, const char *outdir, int blocknum )
   {
   string infname(inputfile);
   int pos = infname.rfind("/");
   infname = infname.substr( pos+1, infname.length() );

   int dot = infname.find(".");
   infname = infname.substr( 0, dot );

   string filename(outdir);
   filename += infname;

   char fnamet[200];
   sprintf( fnamet, "_BL%d_B%d.txt",  m_nBlockLength, blocknum );

   filename +=  string( fnamet );

   int ii, jj, tt, kk, ee;
   int K = m_A.size();
   FILE *fp = fopen( filename.c_str(), "w" );

// 0. Save Header 
fprintf( fp, "********************************\n" );
fprintf( fp, "*	  MUTDP Haplotyper \n" );
fprintf( fp, "********************************\n\n" );

// 4. Save H
fprintf( fp, "[H]\n" );
int I = m_pData->m_numTotalI;
unsigned char **pred_h[2];
pred_h[0] = m_pData->m_Pred_Haplotypes[0];
pred_h[1] = m_pData->m_Pred_Haplotypes[1];
for (ii=0; ii < I; ii++)
{
fprintf( fp, "%d, G%d --- %d, %d	\n", ii, m_pData->m_EthnicGroup[ii], 
m_EqClass[ii][0], m_EqClass[ii][1] );
for (ee=0; ee<2; ee++)
{
for (tt = m_nBlockStart; tt< m_nBlockEnd; tt++)
fprintf( fp, "%d", pred_h[ee][ii][tt] );

fprintf( fp, "\n" );
}
}

fprintf( fp, "\n" );
// 1. Save A
fprintf( fp, "[A]\n" );
fprintf( fp, "ID	Frequency	%%	Haplotype\n" );
for (kk=0; kk < K; kk++)
{
fprintf( fp, "%d	%d	%.6f	", kk, m_nSumClassN[kk], m_nSumClassN[kk]/(float)(2*I) );
for (tt = 0; tt< m_nBlockLength; tt++)
{
fprintf( fp, "%d", m_A[kk][tt] );
}
fprintf( fp, "\n" );
}

fprintf( fp, "\n%2f\n\n", m_gamma );

// 3. Save DP
for (jj=0; jj<m_dp.size(); jj++)
{
fprintf( fp, "[DP-%d]\n", jj );
fprintf( fp, "%2f\n", m_dp[jj].m_alpha );
for (kk=0; kk < K; kk++)
{
fprintf( fp, "%3d ", m_dp[jj].m_NumClassN[kk] );
}
fprintf( fp, "\n" );
}

fclose(fp);
return 1;
}
    */
int MUTDP::Save( const char *inputfile, const char *outdir, int tstart, int tend )
{
    string infname(inputfile);
    int pos = infname.rfind("/");
    infname = infname.substr( pos+1, infname.length() );

    int dot = infname.find(".");
    infname = infname.substr( 0, dot );

    string filename(outdir);
    filename += infname;

    char fnamet[200];
    sprintf( fnamet, "_hap_T%d_%d.txt", tstart, tend );

    filename +=  string( fnamet );

    int ii, tt, kk, ee;
    int K = m_A.size();
    FILE *fp = fopen( filename.c_str(), "w" );

    // 0. Save Header
    fprintf( fp, "*******************************************\n" );
    fprintf( fp, "*  Phased haplotypes from MUTDP Haplotyper \n" );
    fprintf( fp, "*******************************************\n\n" );

    // 4. Save H
    fprintf( fp, "[H]\n" );
    int I = m_pData->m_numTotalI;
    unsigned char **pred_h[2];
    pred_h[0] = m_pData->m_Pred_Haplotypes[0];
    pred_h[1] = m_pData->m_Pred_Haplotypes[1];
    for (ii=0; ii < I; ii++)
    {
	fprintf( fp, "%d, G%d --- %d, %d	\n", ii, m_pData->m_EthnicGroup[ii],
		m_EqClass[ii][0], m_EqClass[ii][1] );
	for (ee=0; ee<2; ee++)
	{
	    for (tt = m_nBlockStart; tt< m_nBlockEnd; tt++)
		fprintf( fp, "%d", pred_h[ee][ii][tt] );

	    fprintf( fp, "\n" );
	}
    }

    fprintf( fp, "\n" );
    // 1. Save A
    fprintf( fp, "[A]\n" );
    fprintf( fp, "ID	Frequency	%%	Haplotype\n" );
    for (kk=0; kk < K; kk++)
    {
	fprintf( fp, "%d	%d	%.6f	", kk, m_NumClassN[kk], m_NumClassN[kk]/(float)(2*I) );
	for (tt = 0; tt< m_nBlockLength; tt++)
	{
	    fprintf( fp, "%d", m_A[kk][tt] );
	}
	fprintf( fp, "\n" );
    }

    fclose(fp);
    return 1;
}


int MUTDP::Sample_Conparam( int numiter_a, int numiter_b  )
{
    int iter, nd, zz, kk;
    double aa, bb, xx;
    double alpha, gamma;

    int K = m_A.size();
    alpha = m_alpha;

    for ( iter = 0 ; iter < numiter_a ; iter++ )
    {
	aa = alpha_a;
	bb = alpha_b;
	//for ( jj = 0 ; jj < J ; jj++ )
	//{
	    //nd = m_pDataIndex.size();
	    nd =  m_pData->m_numTotalI;
	    xx = randbeta( alpha + 1.0, nd );
	    zz = ( drand48() * (alpha + nd) < nd );
	    // numTable_jj:
	    for (kk=0; kk < K; kk++)
	    {
		if ( m_NumClassN[kk] > 0 )
		    aa++;
	    }
	    aa -= zz;
	    bb -= log(xx);
	    //}
	alpha = randgamma(aa) / bb;
    }
    //for ( jj=0; jj < J; jj++)
    //m_dp[jj].m_alpha = alpha;

    gamma = m_gamma;
    for ( iter = 0 ; iter < numiter_b ; iter++ )
    {
	aa = gamma_a;
	bb = gamma_b;

	nd = 0;
	for (kk=0; kk < K; kk++)
	{
	    if ( m_NumClassN[kk] > 0 )
	    {
		aa++;
	    }
	    nd += m_NumClassN[kk];
	}
	xx = randbeta( gamma + 1.0, nd);
	zz = ( drand48() * (gamma + nd) < nd );

	aa -= zz;
	bb -= log(xx);

	gamma = randgamma(aa) / bb;
    }
    m_gamma = gamma;

    return 1;
}

// Update m_A, m_NumClassN from the predicted h
// Start of Ligation process
// and other SS does not have meaning anymore
int MUTDP::GetPredFreq( unsigned char*** h )
{
    int ii, ee, tt, k, jj;

    int		m_bLigationStep = 1;

    m_A.clear();
    //m_nSumClassN.clear();
    m_NumClassN.clear();

    int I = m_pData->m_numTotalI;

    if ( m_EqClass == 0 )
    {
	m_pData->Alloc2DMemory( &m_EqClass, I, 2);
    }


    int		nstart = m_nBlockStart;
    int		nend = m_nBlockEnd;
    int		nlength = nend - nstart;

    for (ii=0; ii<I; ii++)
    {
	for (ee=0; ee<2; ee++)
	{
	    k = Find( &( h[ee][ii][nstart] ), m_A );
	    if ( k < 0 )
	    {
		vector<unsigned char> hnew;
		for (tt=0; tt<nlength; tt++)
		    hnew.push_back( h[ee][ii][nstart+tt] );

		m_A.push_back( hnew );
		m_NumClassN.push_back(0);
		k = m_A.size() - 1;
	    }
	    m_EqClass[ii][ee] = k;
	    //m_nSumClassN[k] += 1;
	    //m_NumClassN[k] += 1;
	    jj = m_pData->m_EthnicGroup[ii];
	    //	to update the SS of each DP
	}
    }

    return 1;
}

//find class k which is equal to H
int MUTDP::Find( unsigned char *h, vector<vector<unsigned char> > &A,
	int exceptK  )
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


int MUTDP::Swap( unsigned char *h1, unsigned char *h2, int nstart, int nend )
{
    int tt;
    unsigned char tmp;
    for (tt = nstart; tt < nend; tt++)
    {
	tmp = h1[tt];
	h1[tt] = h2[tt];
	h2[tt] = tmp;
    }

    return 1;
}

int MUTDP::ResetCumH()
{
    int **cum_h[B][2];
    int ii, bb, tt;

    for (bb=0; bb<B; bb++)
    {
	cum_h[bb][0] = m_pData->m_Cum_Haplotypes[bb][0];
	cum_h[bb][1] = m_pData->m_Cum_Haplotypes[bb][1];
    }

    int I = m_pData->m_numTotalI;
    for (bb=0; bb<B; bb++)
    {
	for (ii=0; ii<I; ii++)
	{
	    for ( tt= m_nBlockStart; tt < m_nBlockEnd; tt++)
	    {
		cum_h[bb][0][ii][tt] = 0;
		cum_h[bb][1][ii][tt] = 0;
	    }
	}
    }

    return 1;
}

int MUTDP::Initialize( haplo2_t h0, int I, int T, int offset, bool cpShiftedRand )
{
    int		ii, tt, ee, it, cc;
    int		nGroups = m_pData->m_nGroups;
    bool	new_class;

    // Geno-Haplo data
    assert( offset + T <= m_pData->m_numT );

    // Geno-Haplo data
    m_nBlockStart = offset;
    m_nBlockEnd = T + offset;
    m_nBlockLength = T;

    int nstart = offset;
    int nend = T + offset;

    int		numBlockT = m_nBlockLength;

    //  constants
    ab_h = beta_h + alpha_h;
    a1 = ( alpha_h ) / (ab_h );
    b1 = ( beta_h ) / (ab_h );
    ab_g = beta_g + alpha_g;
    a2	= (alpha_g) / (ab_g);
    b2	= (beta_g) / (ab_g);
    logB1 = log(B-1);
    logB = log(B);
    log_mu1=log(mu1);
    log_mu2=log(mu2);
    tiny = pow(10, -100); //10^(-100);

    unsigned char ***h = m_pData->m_Haplotypes;
    unsigned char ***g = m_pData->m_Genotypes;
    int	**g_match = m_pData->m_g_match;
    int	**g_miss1 = m_pData->m_g_miss1;
    int	**g_miss2 = m_pData->m_g_miss2;
    int	***h_count = m_pData->m_h_count;

    vector<int> &n = m_NumClassN;
    //vector<int> &sum_mj = m_nSumClassN;

    //////////////////////////
    // clear previous record
    clearSS();
    m_A.clear();
    //////////////////////////

    // initialize with 1 class
    AddClass( 0, numBlockT );

    //govind
    ////////////////////////////////////////////////
    //m_dp.clear();
    // create DP: in case of known ethnic group variable
    //for ( jj = 0; jj < nGroups; jj++ )
    //{
    //	m_dp.push_back( *(new DP()) );
    //}

    // initialize DP ethnic groups
    //for ( jj = 0; jj < nGroups; jj++ )
    //{
    //	m_dp[jj].m_pDataIndex = m_pData->m_DataInGroup[jj];
    //	m_dp[jj].AddClass( 0, numBlockT );
    //}

    if ( m_EqClass != 0 )
    {
	delete [] m_EqClass[0];
	delete [] m_EqClass[1];
	delete [] m_EqClass;
    }
    m_pData->Alloc2DMemory( &m_EqClass, I, 2);

    // initial assignment
    for ( ii = 0; ii < I ; ii++ )
    {
	if ( cpShiftedRand )
	{
	    for ( tt = nstart; tt < nend; tt++)
	    {
		float aa = rand()/(float)RAND_MAX;
		int aa1 = 0;
		if ( aa > 0.5 )		aa1 = 1;

		h[0][ii][tt] = h0[aa1][ii][tt];
		h[1][ii][tt] = h0[1-aa1][ii][tt];
		if ( h[0][ii][tt] == 255 )
		{
		    h[0][ii][tt] = 0;
		    h[1][ii][tt] = 0;
		}
	    }
	}
	else
	{
	    // init h
	    for (tt = nstart; tt < nend; tt++)
	    {
		h[0][ii][tt] = h0[0][ii][tt - nstart];
		h[1][ii][tt] = h0[1][ii][tt - nstart];
		if ( h[0][ii][tt] == 255 )
		{
		    h[0][ii][tt] = 0;
		    h[1][ii][tt] = 0;
		}
	    }
	}

	// init g_match, g_miss1, g_miss2
	for (tt = nstart; tt < nend; tt++)
	{
	    int g_id0 = MIN( g[0][ii][tt], g[1][ii][tt] );
	    int g_id1 = MAX( g[0][ii][tt], g[1][ii][tt] );
	    int h_id0 = MIN( h[0][ii][tt], h[1][ii][tt] );
	    int h_id1 = MAX( h[0][ii][tt], h[1][ii][tt] );

	    if ( g_id0 == h_id0 && g_id1 == h_id1 )
	    {
		g_match[ii][tt] = 1;
		g_miss1[ii][tt] = g_miss2[ii][tt] = 0;
	    }
	    else if ( g_id0 != h_id0 && g_id1 != h_id1 )
	    {
		g_miss2[ii][tt] = 1;
		g_miss1[ii][tt] = g_match[ii][tt] = 0;
	    }
	    else
	    {
		g_miss1[ii][tt] = 1;
		g_miss2[ii][tt] = g_match[ii][tt] = 0;
	    }
	}
    }

    // calculate u[0:2] for each DP
    CalNumClassU( );

    // rand permutation of ii
    vector<int> iivec( I );
    for ( ii = 0; ii < I; ii++ )	iivec[ii] = ii;
    random_shuffle( iivec.begin(), iivec.end() );


    // initialization other variables
    int ** c = m_EqClass;
    for ( it = 0; it < I ; it++ )
    {
	for ( ee = 0; ee < 2; ee++ )
	{

	    ii = iivec[it];
	    int K = m_NumClassN.size();

	    vector<vector<int> >	&l	= m_NumClassL;
	    vector<vector<int> >	*la = m_NumClassLA;
	    int				*u	= m_NumClassU;

	    double alpha0 = 0.7;

	    /// 1. Sample c(i,e)
	    cc = Sample_EqClass_Init( h[ee][ii], l, alpha0);
		//, &FromTopLevel[ii][ee] );

		c[ii][ee] = cc;

	    if ( cc < K )
	    {
		n[ cc ] ++;
		//sum_mj[ cc ] ++; //no need to do it.
		new_class = 0;

	    }
	    else
	    {
		new_class = 1;
		// add class
		AddClass( 1, numBlockT );	// top-level
		K++;
	    }

	    // update LA
	    for ( tt = nstart; tt < nend; tt++)
	    {
		la[ h[ee][ii][tt] ][cc][ tt-nstart ] ++;
	    }

	    // 2. sample ancestor A_cc
	    m_A[cc] = Sample_A( cc, h[ee][ii], new_class );

	    K = m_A.size();

	    for ( tt = nstart; tt < nend; tt++)
	    {
		la[ h[ee][ii][tt] ][cc][ tt-nstart ] --;
	    }

	    // 3. sample individual haplotype h[ee][ii]
	    Sample_H( h[ee][ii], h[1-ee][ii],
		    g[0][ii], g[1][ii], m_A[cc],
		    g_match[ii], g_miss1[ii], g_miss2[ii],
		    n[cc], l[cc], la[0][cc], la[1][cc],
		    h_count[ee][ii],
		    u, I);

	}	// end of e
    }	// end of i

    return 1;
}

int MUTDP::InferHaplotypes( GenoHaploDB* pDB, int nstart, int nend )
{
    // do adaptive number of iterations by checking convergence (default)
    if ( m_bCheckConvg )
    {
	InferHaplotypesAdaptive( pDB, nstart, nend	 );
    }
    // do fixed number of iterations
    else
    {
	Initialize( pDB, nstart, nend );
	Iterate_det_Gibbs_Met( m_nBurninIteration );
	Iterate_cum_Gibbs_Met( m_nCumIteration );
    }

    return 1;
}

int MUTDP::InferHaplotypesAdaptive( GenoHaploDB* pDB, int nstart, int nend )
{
    int ii, tt;

    /// 1. init
    Initialize( pDB, nstart, nend );

    /// 2. burnin iteration
    Iterate_det_Gibbs_Met( m_nBurninIteration );

    int I = m_pData->m_numTotalI;

    int perIter = MIN( 200*m_nThining, m_nCumIteration );
    int Maxiter = MAX( (int) m_nCumIteration / perIter, 1 );
    m_traceDiff.clear();

    /// 3. first cum iteration
    ////   - will do multiple cum iterations until convergence
    Iterate_cum_Gibbs_Met( perIter );

    unsigned char **pred_h[2];
    pred_h[0] = m_pData->m_Pred_Haplotypes[0];
    pred_h[1] = m_pData->m_Pred_Haplotypes[1];
    unsigned char **predh_old[2];
    m_pData->Alloc2DMemory( &(predh_old[0]), I, m_nBlockLength );
    m_pData->Alloc2DMemory( &(predh_old[1]), I, m_nBlockLength );

    for (ii=0; ii<I; ii++)
    {
	memcpy( predh_old[0][ii], &(pred_h[0][ii][nstart]),
		sizeof(unsigned char)*m_nBlockLength );
	memcpy( predh_old[1][ii], &(pred_h[1][ii][nstart]),
		sizeof(unsigned char)*m_nBlockLength );
    }

    int nTotSite = I*m_nBlockLength*2;
    float thr = nTotSite * 0.01;		// set threshold for convergence check
    int nConvBl = 0, iter = 0;
    ////    - repeat cum iteration until convergence
    while ( iter++ < Maxiter )
    {
	Iterate_cum_Gibbs_Met( perIter );
	Sample_Pred();
	int mat1=0, mat2=0;
	for ( ii=0; ii<I; ii++)
	{
	    for ( tt = nstart; tt < nend ; tt++ )
	    {
		if ( pred_h[0][ii][tt] == predh_old[0][ii][tt-nstart] )
		    mat1++;
		if ( pred_h[1][ii][tt] == predh_old[1][ii][tt-nstart] )
		    mat1++;
		if ( pred_h[0][ii][tt] == predh_old[1][ii][tt-nstart] )
		    mat2++;
		if ( pred_h[1][ii][tt] == predh_old[0][ii][tt-nstart] )
		    mat2++;
	    }
	}
	int diff = nTotSite - MAX( mat1, mat2 );

	// converged if below threshold twice
	if ( diff < thr )
	{
	    if ( nConvBl == 0 ) nConvBl++;
	    else break;
	}
	else
	{
	    if ( nConvBl == 1 ) nConvBl = 0;
	}
	// copy old value
	for (ii=0; ii<I; ii++)
	{
	    memcpy( predh_old[0][ii], &(pred_h[0][ii][nstart]),
		    sizeof(unsigned char)*m_nBlockLength );
	    memcpy( predh_old[1][ii], &(pred_h[1][ii][nstart]),
		    sizeof(unsigned char)*m_nBlockLength );
	}
    }

    // compute predictive estimation
    Sample_Pred();

    return 1;
}

int MUTDP::InferHaplotypes(haplo2_t h0, int* popLabel,
	int I, int T, int offset, bool* bDone )
{
    Initialize( h0, I, T, offset );

    float tmpg = m_gamma;
    int tmpDoCP = m_doConparam;
    m_doConparam = 1;

    Iterate_det_Gibbs_Met( m_nBurninIteration, bDone );
    Iterate_cum_Gibbs_Met( m_nCumIteration, bDone );

    m_gamma = tmpg;
    m_doConparam = tmpDoCP;

    return 1;
}

///////////////// estimate mutation rate ////////////////
int MUTDP::EstimateTheta( int iter )
{
    vector<float> theta;
    int K = m_A.size();
    int  kk, tt;

    vector<int> &nn = m_NumClassN;;
    vector<vector<int> > ll = m_NumClassL;
    int T = ll[0].size();

    float meanTheta = m_pData->EstimatedTheta( nn, ll, theta );

    if ( m_pFPTheta )
    {
	fprintf( m_pFPTheta, "%d %.4f	", iter, meanTheta );
	for (tt=0; tt<T; tt++)
	    fprintf( m_pFPTheta, "%.4f	", theta[tt] );

	int numPopK = 0;
	for (kk=0; kk<K; kk++)
	{
	    if ( nn[kk] > 0 )
		numPopK ++;
	}
	fprintf( m_pFPTheta, "%d %d\n", numPopK, K );
    }
    return 1;
}


int MUTDP::LoadData(const char *filename)
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


int MUTDP::DeleteSS()
{
	m_NumClassN.clear();
	m_NumClassL.clear();
	m_NumClassLA[0].clear();
	m_A.clear();
	m_NumClassLA[1].clear();

	return 1;
}

int MUTDP::DeleteClass()
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


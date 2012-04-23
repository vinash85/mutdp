#include <string.h>
#include <cstring>
#include <time.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
//#include <cstdlib.h>

#ifdef _WIN32
#include <direct.h>
#define vbmkdir(x) _mkdir(x)
#else
#define vbmkdir(x) mkdir(x,S_IRWXU|S_IRGRP|S_IXGRP )
#endif 

#include "MUTDP.h"
#include "HapAssembler.h"

using namespace std;

////////////////// for initialization ////////////////////
int		HAP_BL_LENGTH	= 7;
int		HAP_BL_START	= 0;
int		HAP_BL_END		= 10000;
int		HAP_BL_MAX	= 10000;
int		HAP_BL_MIN	= 0;
int		HAP_NUM_BURN_IN = -1; //5000; //50000;
int		HAP_NUM_CUM = -1; // 3000; // 50000;
int		HAP_THINING = -1;
int		HAP_CONPARAM = -1;
int		HAP_CONVCHECK = -1;
const char *HAP_INI_FILE;
char *HAP_OUTDIR = 0;


///////////////////////////////////////////
/// parse argument
///////////////////////////////////////////
void parseArg( int argc, char *argv[] )
{
	if ( argc < 2 )
	{
		printf("[usage]: Haplo.exe inputfile\n");
		exit(-1);
	}

	HAP_INI_FILE = argv[1];
	int counter = 1;

	while (++counter < argc) 
	{
		if ( strcmp( argv[counter], "-bl" ) == 0 ) // block length in Partition-Ligation for long sequences
		{
			counter++;
			HAP_BL_LENGTH = atoi( argv[counter] );			
		} 
		else if ( strcmp ( argv[counter], "-dir" ) == 0 )	// output dir name
        {
              counter++;
              HAP_OUTDIR = argv[counter] ;
        }
		else if ( strcmp ( argv[counter], "-cp" ) == 0 )	// flag for resampling DP hyper parameters
		{
			counter++;
			HAP_CONPARAM = atoi( argv[counter] );			
		}
		else if ( strcmp ( argv[counter], "-conv" ) == 0 )	// flag for checking convergence or doing fixed # of iter
		{
			counter++;
			HAP_CONVCHECK = atoi( argv[counter] );			
		}

		else if ( strcmp ( argv[counter], "-nbr" ) == 0 )	// number of burnin iterations
		{
			counter++;
			HAP_NUM_BURN_IN = atoi( argv[counter] );			
		}
		else if ( strcmp ( argv[counter], "-nc" ) == 0 )	// number of cumulative iterations
		{
			counter++;
			HAP_NUM_CUM = atoi( argv[counter] );		
		}
		else if ( strcmp ( argv[counter], "-nthin" ) == 0 )	// thining interval
		{
			counter++;
			HAP_THINING = atoi( argv[counter] );		
		}
		else if ( strcmp ( argv[counter], "-mbl" ) == 0 )	// max SNP index
		{
			counter++;
			HAP_BL_MAX = atoi( argv[counter] );			
		}
		else if ( strcmp ( argv[counter], "-mnbl" ) == 0 )	// min SNP index
		{
			counter++;
			HAP_BL_MIN = atoi( argv[counter] );			
		}
		else 
		{
			printf( "Unknown option: [%s]\n", argv[counter]);
//		    exit(-1);
		}		
	} // while
}

///////////////////////////////////////////
/// Load ini file
///////////////////////////////////////////
vector<string> LoadIniFile( const char* inifile )
{
	char	tmp[500];
	int		nfile;
	vector<string>	filelist;
	
	FILE *fp = fopen( inifile, "r" );
	if ( fp == 0 )
	{
		printf( "File Not Exist: %s\n", inifile );
		exit(-1);
	}

	fscanf( fp, "%d", &nfile );			// line 1. number of files (populations)
	for (int i=0; i < nfile; i++)		// line 2-n. filename for each population data
	{
		fscanf( fp, "%s", tmp );
		filelist.push_back( tmp );
		printf( "%s\n", filelist[i].c_str());
	}
	
	fclose(fp);
	
	return filelist;
}

///////////////////////////////////////////
/// Make output directory
///////////////////////////////////////////
int MakeDir( const char *dir )
{
	string dname( dir );

	int sep, sti=0;
	string rem(dname);
	string sub;
	while (1)
	{
		sep = rem.find( "/" );
		if ( sep == -1 )
			return 1;

		sti += sep+1;
		sub = dname.substr(0, sti);
		vbmkdir( sub.c_str() );
		rem = rem.substr( sep+1, rem.length() );
	}

	return 1;
}

///////////////////////////////////////////
int main( int argc, char *argv[] )
{

	/////////////////////////////////
	// parse arguments
	parseArg( argc, argv );

	/////////////////////////////////
	// Read ini file
	vector<string>	filelist = LoadIniFile( HAP_INI_FILE );

	//////////////////////////////////
	// Load DB
	GenoHaploDB	DB;
	DB.m_nMaxT = HAP_BL_MAX;
	DB.m_nMinT = HAP_BL_MIN;
	DB.LoadData( filelist );
	HAP_BL_MAX = DB.m_nMaxT;
	int		I = DB.m_numTotalI;
	int		numT = DB.m_numT;
	int		nGroups = filelist.size();

	//////////////////////////////////
	clock_t start, finish;
	start = clock();

	///////////////////////////////////////////////////
	///////// parse file path ////////////
	string outdir( "./output/" );
	string dbname( HAP_INI_FILE );

	int pos = dbname.rfind("/");
	if ( pos > 0 )
			dbname = dbname.substr( pos+1, dbname.length() );
	int dot = dbname.rfind( "." );
	dbname = dbname.substr(0, dot);
		
	if ( HAP_OUTDIR != 0 )
		outdir = HAP_OUTDIR;

	MakeDir( outdir.c_str() );


	////////////////////////////////////////////////
	////	create HDP instance to phase each block
	MUTDP Haplo;
	Haplo.Init( &DB, 0, numT );
		//Burning is transitory starting time of markov chain
	if ( HAP_NUM_BURN_IN > 0 )
		Haplo.m_nBurninIteration = HAP_NUM_BURN_IN;
	if ( HAP_NUM_CUM > 0 )
		Haplo.m_nCumIteration = HAP_NUM_CUM;
	if ( HAP_THINING > 0 )
	    Haplo.m_nThining = HAP_THINING;
	if (  HAP_CONPARAM > 0 )
		Haplo.m_doConparam = HAP_CONPARAM;
	if (  HAP_CONVCHECK > 0 )
		Haplo.m_bCheckConvg = HAP_CONVCHECK;

	////////////////////////////////////////////////
	////	create block-assembler to ligate blockwise result
	HapAssembler Assembler( outdir, dbname ) ;

	Assembler.Initialize( &DB, &Haplo );
	vector<haplo2_t> &arrData = Assembler.m_arrData;
	vector<haplo2_t> &arrPred = Assembler.m_arrPairwisePred;


	int res = numT % HAP_BL_LENGTH;
	int	numBlocks = (numT - res) / HAP_BL_LENGTH;
	if ( res !=0 ) 
			numBlocks++;

	
	int		nstart, nend, nlength;
	///////////////////////////////////////
	///		1. Block-wise inference  using HDP haplotyper (Partition)
	///////////////////////////////////////
	//TODO change this numBlcoks = 1
	for ( int nb = 0; nb < numBlocks; nb++ )
	{
		nend = MIN( ( nb + 1 ) * HAP_BL_LENGTH, numT );
		nstart = MAX( nend - HAP_BL_LENGTH, 0);
		nlength = nend - nstart;
		if ( HAP_BL_LENGTH != nlength )
		{	
			HAP_BL_LENGTH = nlength;
			cout << "block length not equal to " << HAP_BL_LENGTH << endl;
		}	

		Haplo.InferHaplotypes( &DB, nstart, nend );

		// save intermediate result ( just for checking, can be commented out )
		Haplo.Save( dbname.c_str(), outdir.c_str(), nstart, nend );

		// add the blockwise result to the output array
		haplo2_t	*phh = new haplo2_t;
		DB.AllocHaplo2( phh, I, nlength );
		DB.CopyHaplo( nstart, nend, phh );
		arrData.push_back( *phh );
	}
	
	///////////////////////////////////////////////
	///////// 2. Ligation 
	///////////////////////////////////////////////
	vector<int> T1( numBlocks, HAP_BL_LENGTH );		// length of each block
	vector<int> nOverlap( numBlocks, 0 );			// overlap between adjacent block
	if ( res != 0 )
		nOverlap[ numBlocks - 2 ] = HAP_BL_LENGTH - res;

	Assembler.LigateHaplotypes( arrData, T1, nOverlap, arrPred );


	//////////// save final result ///////////////
	// output filename
	string ofilename( outdir );
	ofilename += dbname + "_final.txt";
	Assembler.SaveResult( ofilename.c_str() );
	
	finish = clock();
	
	printf("Elapsed time: %.2f min\n", (double)(finish-start)/CLOCKS_PER_SEC/60.0 );

	return 1;
}

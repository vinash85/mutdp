#include <vector>
#include <math.h>
#include <stdlib.h>
#include <cstdlib>
#include <iostream>



using namespace std;

// normalize vector to have unit length
int normalize( vector<double> &vec )
{
	int n = vec.size();
	double norm = 0;
	vector<double>::iterator it;
	
	double max = -INF;
	// get L1 norm
	for (it = vec.begin(); it!=vec.end(); ++it)
	{
		double its = (*it);
		norm += its;
		if ( its > max )
			max = its;
	}
//	norm = sqrt( norm );
	
	// normalize
	for (it = vec.begin(); it!=vec.end(); ++it)
	{
		(*it) = (*it)/norm;
	}
	
	return 1;
}


// normalize vector to have unit length
template <class T>
T	sum( vector<T> &invec )
{
	T sumval = 0;
	int nn = invec.size();
	for (int it = 0; it < nn; it++)
	{
		sumval += invec[it];
	}
	return sumval;
}

template <class T>
int sample_discrete( vector<T> &pred )
{
	int c=0;
	int n = pred.size();
	double R = rand()/(double)RAND_MAX;
	
	int ind = 0;
	double cumprob = 0;
//	for (vector<T>::iterator it = pred.begin(); it!=pred.end(); ++it, ++ind )
	for (int kk=0; kk<pred.size(); kk++) 
	{
		cumprob += pred[kk]; //(*it);
		if ( R <= cumprob ) 
		{
			c = kk;
			break;
		}		
	}
	return c;
}

// lgamma.cpp -- log gamma function of real argument.
//      Algorithms and coefficient values from "Computation of Special
//      Functions", Zhang and Jin, John Wiley and Sons, 1996.
//
//  (C) 2003, C. Bond. All rights reserved.
//
//  Returns log(gamma) of real argument.
//  NOTE: Returns 1e308 if argument is 0 or negative.
//
// #include <math.h>
double lgamma(double x)
{
    double x0,x2,xp,gl,gl0;
    int n,k;
    static double a[] = {
        8.333333333333333e-02,
       -2.777777777777778e-03,
        7.936507936507937e-04,
       -5.952380952380952e-04,
        8.417508417508418e-04,
       -1.917526917526918e-03,
        6.410256410256410e-03,
       -2.955065359477124e-02,
        1.796443723688307e-01,
       -1.39243221690590};
    
    x0 = x;
    if (x <= 0.0) return 1e308;
    else if ((x == 1.0) || (x == 2.0)) return 0.0;
    else if (x <= 7.0) {
        n = (int)(7-x);
        x0 = x+n;
    }
    x2 = 1.0/(x0*x0);
    xp = 2.0*3.141592 ;
    gl0 = a[9];
    for (k=8;k>=0;k--) {
        gl0 = gl0*x2 + a[k];
    }
    gl = gl0/x0+0.5*log(xp)+(x0-0.5)*log(x0)-x0;
    if (x <= 7.0) {
        for (k=1;k<=n;k++) {
            gl -= log(x0-1.0);
            x0 -= 1.0;
        }
    }
    return gl;
}
        

///////////////////////////////////////////////
// From Y-W 
#ifndef RANDUTILS
#define RANDUTILS

double drand48()
{
    double randnum = rand();
    return randnum/(double)RAND_MAX;   
}

double randgamma(double rr) {
  double bb, cc, dd;
  double uu, vv, ww, xx, yy, zz;

  if ( rr <= 0.0 ) {
    /* Not well defined, set to zero and skip. */
    return 0.0;
  } else if ( rr == 1.0 ) {
    /* Exponential */
    return - log(drand48());
  } else if ( rr < 1.0 ) {
    /* Use Johnks generator */
    cc = 1.0 / rr;
    dd = 1.0 / (1.0-rr);
    while (1) {
      xx = pow(drand48(), cc);
      yy = xx + pow(drand48(), dd);
      if ( yy <= 1.0 ) {
		  double dr = drand48();
		  if ( dr == 0 )
			  dr = 0.000001;
		
		  return -log( dr ) * xx / yy;
      }
    }
  } else { /* rr > 1.0 */
    /* Use bests algorithm */
    bb = rr - 1.0;
    cc = 3.0 * rr - 0.75;
    while (1) {
      uu = drand48();
      vv = drand48();
      ww = uu * (1.0 - uu);
      yy = sqrt(cc / ww) * (uu - 0.5);
      xx = bb + yy;
      if (xx >= 0) {
        zz = 64.0 * ww * ww * ww * vv * vv;
        if ( ( zz <= (1.0 - 2.0 * yy * yy / xx) ) ||
             ( log(zz) <= 2.0 * (bb * log(xx / bb) - yy) ) ) {
          return xx;
        }
      }
    }
  }
}

int randnumtable(double alpha, int numdata) {
  int ii, numtable;

  if ( numdata == 0 ) {
    return 0;
  } else {
    numtable = 1;
    for ( ii = 1 ; ii < numdata ; ii++ ) {
      if ( drand48() < alpha / (ii+alpha) ) numtable++;
    }
    return numtable;
  }
} 

void randdir(double *pi, double *alpha, int veclength, int skip) {
  double *pi2, *piend;
  double sum;

  sum = 0.0;
  piend = pi + veclength*skip;
  for ( pi2 = pi ; pi2 < piend ; pi2 += skip) {
    sum += *pi2 = randgamma(*alpha);
    alpha += skip;
  }
  for ( pi2 = pi ; pi2 < piend ; pi2 += skip) {
    *pi2 /= sum;
  }
}

double randbeta(double aa, double bb) {
  aa = randgamma(aa);
  bb = randgamma(bb);
  return aa/(aa+bb);
}

int randmult(double *pi, int veclength, int skip) {
  double *pi2, *piend;
  double sum = 0.0, mass;
  int cc = 0;

  piend = pi + veclength*skip;
  for ( pi2 = pi ; pi2 < piend ; pi2 += skip )
    sum += *pi2;
  mass = drand48() * sum;
  while (1) {
    mass -= *pi;
    if ( mass <= 0.0 ) break;
    pi += skip;
    cc ++;
  }
  return cc;
}

int randuniform(int numvalue) {
 	int val = (int)floor( drand48() * numvalue ); 
	return val ;
}

#endif

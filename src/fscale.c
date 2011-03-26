/**********************************************************************
 * 
 * fscale.c
 *
 * copyright (c) 2011, Karl W Broman
 *
 * last modified Mar, 2011
 * first written Mar, 2011
 *
 *     This program is free software; you can redistribute it and/or
 *     modify it under the terms of the GNU General Public License,
 *     version 3, as published by the Free Software Foundation.
 * 
 *     This program is distributed in the hope that it will be useful,
 *     but without any warranty; without even the implied warranty of
 *     merchantability or fitness for a particular purpose.  See the GNU
 *     General Public License, version 3, for more details.
 * 
 *     A copy of the GNU General Public License, version 3, is available
 *     at http://www.r-project.org/Licenses/GPL-3
 *
 * C functions for the R/lineup package
 *
 * Contains: R_fscale, fscale
 *
 **********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include "fscale.h"
#include "util.h"

/* for calling fscale() from R */
void R_fscale(int *nrow, int *ncol, double *x)
{
  double **X;

  reorg_dmatrix(*nrow, *ncol, x, &X);
  
  fscale(*nrow, *ncol, X);
}

/**********************************************************************
 * 
 * fscale: standardize columns of a matrix 
 *
 * Notes:  X indexed as X[col][row] 
 *
 *         The one-pass method can have a lot of round-off error,  
 *         but it is quick.
 *
 **********************************************************************/
void fscale(int nrow, int ncol, double **X)
{
  int i, j, n;
  double sum, sumsq, first, diff;
  
  for(j=0; j<ncol; j++) {
    sum = sumsq = 0.0;
    n = 0;
    first=NA_REAL;
    for(i=0; i<nrow; i++) {
      if(R_FINITE(X[j][i])) {
	n++;
	if(!R_FINITE(first)) first = X[j][i]; /* first non-missing value */
	else {
	  /* sum(x) and sum(x*x) with x centered at first non-missing value*/
	  sum += (diff=(X[j][i]-first)); 
	  sumsq += (diff*diff);
	}
      }
    }
    if(n > 1) { /* if n < 2, do nothing */
      sumsq = sqrt((sumsq - (sum*sum)/(double)n)/(double)(n-1));
      sum /= (double)n;
      for(i=0; i<nrow; i++) 
	if(R_FINITE(X[j][i])) 
	  X[j][i] = (X[j][i] - sum - first)/(sumsq);
    }
  }
}

/* end of fscale.c */

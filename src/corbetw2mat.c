/**********************************************************************
 * 
 * corbewt2mat.c
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
 * Contains: R_corbetw2mat_paired, corbetw2mat_paired,
 *           R_corbetw2mat_unpaired_lr, corbetw2mat_unpaired_lr
 *           R_corbetw2mat_unpaired_best, corbetw2mat_unpaired_best
 *
 **********************************************************************/

#include <R.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "corbetw2mat.h"
#include "util.h"
#include "mystandardize.h"

void R_corbetw2mat_paired(int *nrow, int *ncol, double *x, double *y,
			  int *scaled, double *cor)
{
  double **X, **Y;

  reorg_dmatrix(*nrow, *ncol, x, &X);
  reorg_dmatrix(*nrow, *ncol, y, &Y);
  
  corbetw2mat_paired(*nrow, *ncol, X, Y, *scaled, cor);
}

void corbetw2mat_paired(int nrow, int ncol, double **X, double **Y,
			int scaled, double *cor)
{
  int i, j, n;
  double temp;

  if(!scaled) { /* scale columns to have mean 0 and SD 1 */
    mystandardize(nrow, ncol, X);
    mystandardize(nrow, ncol, Y);
  } 
    
  for(j=0; j<ncol; j++) {
    temp=0.0; 
    n=0;
    for(i=0; i<nrow; i++) {
      if(R_FINITE(X[j][i]) && R_FINITE(Y[j][i])) {
	temp += (X[j][i] * Y[j][i]);
	n++;
      }
    }
    if(n > 1) cor[j] = temp/(double)(n-1);
  }

}


void R_corbetw2mat_unpaired_lr(int *nrow, int *ncolx, double *x, 
			       int *ncoly, double *y,
			       int *scaled, double *cor, 
			       int *index)
{
  double **X, **Y;

  reorg_dmatrix(*nrow, *ncolx, x, &X);
  reorg_dmatrix(*nrow, *ncoly, y, &Y);
  
  corbetw2mat_unpaired_lr(*nrow, *ncolx, X, *ncoly, Y, *scaled, cor, index);
}

void corbetw2mat_unpaired_lr(int nrow, int ncolx, double **X, 
			     int ncoly, double **Y,
			     int scaled, double *cor, int *index)
{
  int i, jx, jy, n, theindex;
  double temp, themax;

  if(!scaled) { /* scale columns to have mean 0 and SD 1 */
    mystandardize(nrow, ncolx, X);
    mystandardize(nrow, ncoly, Y);
  } 
    
  for(jx=0; jx<ncolx; jx++) {
    theindex = NA_INTEGER;
    themax = -2.0;
    for(jy=0; jy<ncoly; jy++) {
      temp=0.0; 
      n=0;
    
      for(i=0; i<nrow; i++) {
	if(R_FINITE(X[jx][i]) && R_FINITE(Y[jy][i])) {
	  temp += (X[jx][i] * Y[jy][i]);
	  n++;
	}
      }
      if(n > 1) {
	temp /= (double)(n-1);
	if(temp > themax) {
	  themax = temp;
	  theindex = jy;
	}
      }

    } /* end loop over col of y */	  
    cor[jx] = themax;
    index[jx] = theindex+1;

  } /* end loop over col of x */
}



void R_corbetw2mat_unpaired_best(int *nrow, int *ncolx, double *x, 
				 int *ncoly, double *y,
				 int *scaled, double *cor, 
				 int *xindex, int *yindex, 
				 int *numpairs, double *corthresh)
{
  double **X, **Y;

  reorg_dmatrix(*nrow, *ncolx, x, &X);
  reorg_dmatrix(*nrow, *ncoly, y, &Y);
  
  corbetw2mat_unpaired_best(*nrow, *ncolx, X, *ncoly, Y, *scaled, cor, 
			    xindex, yindex, numpairs, *corthresh);
}

void corbetw2mat_unpaired_best(int nrow, int ncolx, double **X, 
			       int ncoly, double **Y,
			       int scaled, double *cor, 
			       int *xindex, int *yindex, 
			       int *numpairs, double corthresh)
{
  int i, jx, jy, n;
  double temp;

  if(!scaled) { /* scale columns to have mean 0 and SD 1 */
    mystandardize(nrow, ncolx, X);
    mystandardize(nrow, ncoly, Y);
  } 
    
  *numpairs = 0;
  for(jx=0; jx<ncolx; jx++) {
    for(jy=0; jy<ncoly; jy++) {
      temp=0.0; 
      n=0;
    
      for(i=0; i<nrow; i++) {
	if(R_FINITE(X[jx][i]) && R_FINITE(Y[jy][i])) {
	  temp += (X[jx][i] * Y[jy][i]);
	  n++;
	}
      }
      if(n > 1) {
	temp /= (double)(n-1);
	if(temp > corthresh) {
	  cor[*numpairs] = temp;
	  xindex[*numpairs] = jx+1;
	  yindex[*numpairs] = jy+1;
	  (*numpairs)++;
	}
      }

    } /* end loop over col of y */	  
  } /* end loop over col of x */
}


/* end of corbetw2mat.c */

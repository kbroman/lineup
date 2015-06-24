/**********************************************************************
 *
 * rmsd.c
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
 * Contains: R_rmsd, rmsd
 *
 **********************************************************************/

#include <R.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "rmsd.h"
#include "util.h"

void R_rmsd(int *nrow, int *ncolx, double *x,
            int *ncoly, double *y, double *d,
            int *symmetric)
{
    double **X, **Y, **D;

    reorg_dmatrix(*nrow, *ncolx, x, &X);
    reorg_dmatrix(*nrow, *ncoly, y, &Y);
    reorg_dmatrix(*ncolx, *ncoly, d, &D);

    rmsd(*nrow, *ncolx, X, *ncoly, Y, D, *symmetric);
}

void rmsd(int nrow, int ncolx, double **X, int ncoly, double **Y,
          double **D, int symmetric)
{
    int i, j, k, n;
    double temp;

    if(symmetric) {
        for(i=0; i<ncolx-1; i++) {
            for(j=(i+1); j<ncoly; j++) {
                D[i][j] = 0.0;
                n = 0;
                for(k=0; k<nrow; k++) {
                    if(R_FINITE(X[i][k]) && R_FINITE(Y[j][k])) {
                        temp = (X[i][k] - Y[j][k]);
                        D[i][j] += temp*temp;
                        n++;
                    }
                }
                D[j][i] = D[i][j] = sqrt(D[i][j]/(double)n);
            }
        }
    }

    else {
        for(j=0; j<ncoly; j++) {
            for(i=0; i<ncolx; i++) {
                D[j][i] = 0.0;
                n = 0;
                for(k=0; k<nrow; k++) {
                    if(R_FINITE(X[i][k]) && R_FINITE(Y[j][k])) {
                        temp = (X[i][k] - Y[j][k]);
                        D[j][i] += temp*temp;
                        n++;
                    }
                }
                D[j][i] = sqrt(D[j][i]/(double)n);
            }
        }
    }

}

/* end of rmsd.c */

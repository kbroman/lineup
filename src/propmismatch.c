/**********************************************************************
 *
 * propmismatch.c
 *
 * copyright (c) 2011, Karl W Broman
 *
 * last modified Apr, 2011
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
 * Contains: R_propmismatch, propmismatch
 *
 **********************************************************************/

#include <R.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "propmismatch.h"
#include "util.h"

void R_propmismatch(int *nrow, int *ncolx, int *x,
                    int *ncoly, int *y, double *wts,
                    double *prop, double *denom)
{
    double **Prop, **Denom;
    int **X, **Y;

    reorg_imatrix(*nrow, *ncolx, x, &X);
    reorg_imatrix(*nrow, *ncoly, y, &Y);
    reorg_dmatrix(*ncolx, *ncoly, prop, &Prop);
    reorg_dmatrix(*ncolx, *ncoly, denom, &Denom);

    propmismatch(*nrow, *ncolx, X, *ncoly, Y, wts, Prop, Denom);
}

void propmismatch(int nrow, int ncolx, int **X, int ncoly, int **Y,
                  double *wts, double **Prop, double **Denom)
{
    int i, j, k;
    double temp1, temp2;

    for(j=0; j<ncoly; j++) {
        for(i=0; i<ncolx; i++) {
            temp1 = temp2 = 0.0;
            for(k=0; k<nrow; k++) {

                /* have a bit of trouble regarding NAs being converted to smallest integer */
                if(R_FINITE(X[i][k]) && X[i][k]>INT_MIN &&
                   R_FINITE(Y[j][k]) && Y[j][k]>INT_MIN) {
                    temp2 += wts[k];
                    temp1 += ((double)(X[i][k] != Y[j][k]) * wts[k]);
                }
            }
            Denom[j][i] = temp2;
            if(temp2 > 0) Prop[j][i] = temp1/temp2;
        }
    }
}

/* end of propmismatch.c */

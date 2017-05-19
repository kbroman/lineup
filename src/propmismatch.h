/**********************************************************************
 *
 * propmismatch.h
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
 *     at https://www.r-project.org/Licenses/GPL-3
 *
 * C functions for the R/lineup package
 *
 * Contains: R_propmismatch, propmismatch
 *
 **********************************************************************/

void R_propmismatch(int *nrow, int *ncolx, int *x,
                    int *ncoly, int *y, double *wts,
                    double *prop, double *denom);

void propmismatch(int nrow, int ncolx, int **X, int ncoly, int **Y,
                  double *wts, double **Prop, double **Denom);

/* end of propmismatch.h */

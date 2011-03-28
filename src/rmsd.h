/**********************************************************************
 * 
 * rmsd.h
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

void R_rmsd(int *nrow, int *ncolx, double *x, 
		int *ncoly, double *y, double *d,
		int *symmetric);

void rmsd(int nrow, int ncolx, double **X, int ncoly, double **Y,
	      double **D, int symmetric);

/* end of rmsd.h */

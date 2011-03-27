/**********************************************************************
 * 
 * corbewt2mat.h
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
 *           R_corbetw2mat_unpaired_all, corbetw2mat_unpaired_all
 *           R_corbetw2mat_self, corbetw2mat_self
 *
 **********************************************************************/


void R_corbetw2mat_paired(int *nrow, int *ncol, double *x, double *y,
			  int *scaled, double *cor);

void corbetw2mat_paired(int nrow, int ncol, double **X, double **Y,
			int scaled, double *cor);


void R_corbetw2mat_unpaired_lr(int *nrow, int *ncolx, double *x, 
			       int *ncoly, double *y,
			       int *scaled, double *cor, 
			       int *index);

void corbetw2mat_unpaired_lr(int nrow, int ncolx, double **X, 
			     int ncoly, double **Y,
			     int scaled, double *cor, int *index);

void R_corbetw2mat_unpaired_best(int *nrow, int *ncolx, double *x, 
				 int *ncoly, double *y,
				 int *scaled, double *cor, 
				 int *xindex, int *yindex, 
				 int *numpairs, double *corthresh);

void corbetw2mat_unpaired_best(int nrow, int ncolx, double **X, 
			       int ncoly, double **Y,
			       int scaled, double *cor, 
			       int *xindex, int *yindex, 
			       int *numpairs, double corthresh);

void R_corbetw2mat_unpaired_all(int *nrow, int *ncolx, double *x, 
				int *ncoly, double *y,
				int *scaled, double *cor);

void corbetw2mat_unpaired_all(int nrow, int ncolx, double **X, 
			      int ncoly, double **Y,
			      int scaled, double **Cor);

void R_corbetw2mat_self(int *nrow, int *ncol, double *x, 
			int *scaled, double *cor);

void corbetw2mat_self(int nrow, int ncol, double **X, 
		      int scaled, double **Cor);

/* end of corbetw2mat.h */

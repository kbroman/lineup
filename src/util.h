/**********************************************************************
 *
 * util.h
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
 * These are utility functions
 *
 * Contains: reorg_dmatrix, reorg_imatrix
 *
 **********************************************************************/

/**********************************************************************
 *
 * reorg_dmatrix
 *
 * Reorganize a singly-indexed matrix of doubles so it is a
 * doubly-indexed array rather than a single long vector
 *
 * Afterwards, X indexed like X[col][row]
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/
void reorg_dmatrix(int nrow, int ncol, double *x, double ***X);

/**********************************************************************
 *
 * reorg_imatrix
 *
 * Reorganize a singly-indexed matrix of ints so it is a
 * doubly-indexed array rather than a single long vector
 *
 * Afterwards, X indexed like X[col][row]
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/
void reorg_imatrix(int nrow, int ncol, int *x, int ***X);

/* end of util.h */

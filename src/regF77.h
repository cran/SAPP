/*
*
*  SAPP : Statistical Analysis of Point Processes
*  Copyright (C) 2010    The Institute of Statistical Mathematics
*
*    This program is free software; you can redistribute it and/or modify
*    it under the terms of the GNU General Public License as published by
*    the Free Software Foundation; either version 2 of the License, or
*    (at your option) any later version.
*
*    This program is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU General Public License for more details.
*
*    You should have received a copy of the GNU General Public License
*    along with this program; if not, write to the Free Software
*    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA.
*
*  ismrp at grp.ism.ac.jp
*/

#include <R.h>
#include <Rinternals.h>
#include <libintl.h>

#define _(String) (String)

/* Fortran : */

void F77_NAME(eptrenf)(double *y, double *t, int *n, int *nfunct, int *npara,
              int *nsub, double *cycle, double *xa, double *aic, double *aicm,
              int *morder, double *xval, double *fval, double *x, double *g,
              int *id, double *rmd, double *eee, int *nl, int *nmax, int *np,
              int *nlmax);

void F77_NAME(etarppf)(double *time, double *mag, double *refer, int *nn,
              double *param, int *np, double *zts, double *zte, double *tstart,
              double *x, int *ntst);

void F77_NAME(etasapf)(double *xx, double *mag, int *nn, double *refer,
              double *thresh, double *param, int *np, double *zts, double *zte,
              double *tstart, int *nfunct, int *app, double *f, double *x,
              double *g, double *aic2, int *id, double *ee, double *x1, int *nl,
              int *nlmax);

void F77_NAME(etasimf)(int *ic, double *bvalue, double *tstart, int *nd,
              double *thresh, double *refer, double *a, double *b, double *c,
              double *d, double *p, double *mag1, double *time1, double *mag2,
              double *time2, double *probx);

void F77_NAME(linlinf)(int *n, double *x, int *opt, double *t, int *nn, int *mm,
              double *xx, double *yy, int *kkx, int *kky, int *kmax, int *kkc,
              int *kkt,  int *nlmax, double *x1, double *x2, double *aic,
              double *f, double *prb, double *r1, double *rwx, double *rwy,
              double *phs, double *px, double *pg, int *id, double *rmd,
              double *eee, int *nl, int *ier);

void F77_NAME(linsimf)(int *kxx, int *kxy, int *kxz, double *t, double *c,
              double *d, double *ax, double *ay, double *at, double *yy,
              int *mm, double *ptmax, int *kmax, double *xx, int *ii1, int *jj1,
              double *err, int *ier);

void F77_NAME(momorif)(double *y, int *n, double *pai, int *np, double *zts,
              double *zte, int *ncount, int *nfunct, double *ff, double *x,
              double *g, double *pa, double *ahaic, double *t0, double *ti,
              double *ak, double *c, double *p, double *cls, int *id,
              double *rmd, double *x1, double *h, double *hf, int *nl,
              int *nlmax);

void F77_NAME(pgraphf)(int *nfunct, int *isi, double *zd, int *nn, int *npoint,
              double *days, double *h, double *delta, double *dmax, int *kmax,
              double *xtau, double *y, int *kn, double *xl, double *xx,
              double *ydev, double *ui, double *ncn, double *sui, double *xp,
              double *xrate, double *dlt, double *xtime, double *yvar,
              double *sigma, int *k, int *ier);

void F77_NAME(ptspecf)(double *data, int *n, double *intval, double *pprd,
              double *prdm, double *prd, int *nfre, int *nt, int *is,
              double *prb, double *r1, double *rwx, double *rwy, double *phs,
              double *wt, double *ht, double *w, double *h, double *g);

void F77_NAME(respoif)(double *time, double *mag, double *dep, double *xp,
              double *yp, int *nd, double *param, int *np, double *zts,
              double *zte, double *tstart, double *thresh, double *mag1,
              double *dep1, double *xp1, double *yp1, int *ntstar, double *xx,
              double *x, int *nn);

void F77_NAME(simbvhf)(int *kxx, int *kxy, int *kxz, int *kyx, int *kyy,
              int *kyz, double *t, double *coxx, double *coxy, double *coyx,
              double *coyy, double *axx, double *axy, double *axz, double *ayx,
              double *ayy, double *ayz, double *ptxmax, double *ptymax,
              int *kmax, double *x, double *y, int *ii1, int *jj1, double *err,
              int *ier, int *nnmax, int *mmmax);

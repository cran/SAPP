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

#include "regF77.h"
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

/* .Fortran calls */

static const R_FortranMethodDef FortEntries[] = {
    {"eptrenf",    (DL_FUNC) &F77_NAME(eptrenf),    22},
	{"etarppf",    (DL_FUNC) &F77_NAME(etarppf),    11},
	{"etasapf",    (DL_FUNC) &F77_NAME(etasapf),    21},
    {"etasimf",    (DL_FUNC) &F77_NAME(etasimf),    16},
    {"linlinf",    (DL_FUNC) &F77_NAME(linlinf),    30},
    {"linsimf",    (DL_FUNC) &F77_NAME(linsimf),    18},
    {"momorif",    (DL_FUNC) &F77_NAME(momorif),    26},
    {"pgraphf",    (DL_FUNC) &F77_NAME(pgraphf),    27},
    {"ptspecf",    (DL_FUNC) &F77_NAME(ptspecf),    19},
    {"respoif",    (DL_FUNC) &F77_NAME(respoif),    20},
    {"simbvhf",    (DL_FUNC) &F77_NAME(simbvhf),    28},
    {NULL, NULL, 0}
};

void attribute_visible R_init_SAPP(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

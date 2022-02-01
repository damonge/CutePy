#ifndef _CUTE_
#define _CUTE_

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <complex.h>
#include <omp.h>
#include <fitsio.h>
#include <chealpix.h>

#define MIN(a, b) (((a) < (b)) ? (a) : (b)) //Minimum of two numbers

void cute_correlation_scaled(long ngal,double *pos,long nside,double *fld,double *msk,
			     double thmin,double thmax,int nth,int do_log,
			     double *hf_th,double *hm_th);
void cute_correlation_scaled_2D(long ngal,double *pos,long nside,double *fld,double *msk,
				double thmin,double thmax,int nth,int na,
				int do_log,double *hf_th,double *hm_th);

void cute_line_correlation(long nside, double **fld, double *msk,
			   double xmin, double xmax, int nx, int do_log,
			   double *hf_th, double *hm_th, int per_bin);
#endif //_CUTE_

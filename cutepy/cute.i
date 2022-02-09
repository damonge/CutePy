%module cutelib

%{
#define SWIG_FILE_WITH_INIT
#include "../src/cute.h"
%}

%include "numpy.i"
%init %{
  import_array();
%}

%include "../src/cute.h"

%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* dout, int ndout)};
%apply (int DIM1,double *IN_ARRAY1) {(int npos, double *pos),
                                     (int npix, double *fld),
                                     (int npixs, double *flds),
                                     (int npix2, double *msk)};

%inline %{
void cute_scaled_correlation_wrap(int npos, double *pos,
				  int npix, double *fld,
				  int npix2, double *msk,
				  double thmax, int nth, int nside,
				  double *dout, int ndout)
{
  int ngal=npos/4;
  double *hf_th = dout;
  double *hm_th = &(dout[nth]);
  
  cute_correlation_scaled(ngal, pos, (long)nside, fld, msk,
			  0., thmax, nth, 0,
			  hf_th, hm_th);
}

void cute_scaled_correlation_2D_wrap(int npos, double *pos,
				     int npix, double *fld,
				     int npix2, double *msk,
				     double thmax, int nth, int na,
				     int nside,
				     double *dout, int ndout)
{
  int ngal=npos/4;
  double *hf_th = dout;
  double *hm_th = &(dout[nth*na]);
  
  cute_correlation_scaled_2D(ngal, pos, (long)nside, fld, msk,
			     0., thmax, nth, na, 0,
			     hf_th, hm_th);
}

void cute_line_correlation_wrap(int npixs, double *flds,
				int npix2, double *msk,
				double thmax, int nth, int nside,
				int per_bin, double *dout, int ndout)
{
  double **fld;
  int i, nfld=1;
  double *hf_th = dout;
  double *hm_th = &(dout[nth]);

  if(per_bin)
    nfld = nth;

  fld = malloc(nfld*sizeof(double *));
  for(i=0;i<nfld;i++)
    fld[i] = &(flds[i*npix2]);

  cute_line_correlation((long)nside, fld, msk,
			0., thmax, nth, 0, hf_th,
			hm_th, per_bin);
  free(fld);
 }

void cute_correlation_wrap(int npix, double *fld,
			   int npix2, double *msk,
			   double thmax, int nth, int nside,
			   double *dout, int ndout)
{
  double *hf_th = dout;
  double *hm_th = &(dout[nth]);

  cute_correlation((long)nside, fld, msk,
		   0., thmax, nth, 0,
		   hf_th, hm_th);
 }
%}

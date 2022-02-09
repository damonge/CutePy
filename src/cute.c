#include "cute.h"

static void report_error(int level,char *fmt,...)
{
  va_list args;
  char msg[256];

  va_start(args,fmt);
  vsprintf(msg,fmt,args);
  va_end(args);
  
  if(level) {
    fprintf(stderr," Fatal error: %s",msg);
    exit(level);
  }
  else
    fprintf(stderr," Warning: %s",msg);
}

static void *my_malloc(size_t size)
{
  void *outptr=malloc(size);
  if(outptr==NULL) report_error(1,"Out of memory\n");

  return outptr;
}

static void *my_calloc(size_t nmemb,size_t size)
{
  void *outptr=calloc(nmemb,size);
  if(outptr==NULL)
    report_error(1,"Out of memory\n");

  return outptr;
}

static long he_nside2npix(long nside)
{
  return 12*nside*nside;
}

static void he_pix2vec_ring(long nside, long ipix, double *vec)
{
  pix2vec_ring(nside,ipix,vec);
}

static int nint_he(double n)
{
  if(n>0) return (int)(n+0.5);
  else if(n<0) return (int)(n-0.5);
  else return 0;
}

static int modu_he(int a,int n)
{
  int moda=a%n;
  if(moda<0)
    moda+=n;
  return moda;
}

static void he_in_ring(int nside,int iz,double phi0,double dphi,
		       int *listir,int *nir)
{
  int take_all,conservative,to_top;
  int npix;
  int ncap;
  int diff,jj;
  int nr,nir1,nir2,ir,kshift;
  int ipix1,ipix2;
  int ip_low,ip_hi;
  int nir_here;
  double phi_low,phi_hi,shift;

  conservative=1;//Do we take intersected pixels which
                 //centers do not fall within range?
  take_all=0;//Take all pixels in ring?
  to_top=0;
  npix=(12*nside)*nside;
  ncap=2*nside*(nside-1); //#pixels in north cap
  nir_here=*nir;
  *nir=0;

  phi_low=phi0-dphi-(int)((phi0-dphi)/(2*M_PI))*2*M_PI;
  phi_hi=phi0+dphi-(int)((phi0+dphi)/(2*M_PI))*2*M_PI;
  if(fabs(dphi-M_PI)<1E-6) take_all=1;

  //Identifies ring number
  if((iz>=nside)&&(iz<=3*nside)) {//equatorial
    ir=iz-nside+1;
    ipix1=ncap+4*nside*(ir-1); //Lowest pixel number
    ipix2=ipix1+4*nside-1; //Highest pixel number
    kshift=modu_he(ir,2);
    nr=4*nside;
  }
  else {
    if(iz<nside) {//North pole
      ir=iz;
      ipix1=2*ir*(ir-1);
      ipix2=ipix1+4*ir-1;
    }
    else {//South pole
      ir=4*nside-iz;
      ipix1=npix-2*ir*(ir+1);
      ipix2=ipix1+4*ir-1;
    }
    nr=4*ir;
    kshift=1;
  }

  //Constructs the pixel list
  if(take_all) {
    *nir=ipix2-ipix1+1;
    if(*nir>nir_here)
      report_error(1,"Not enough memory in listir\n");
    for(jj=0;jj<(*nir);jj++)
      listir[jj]=ipix1+jj;

    return;
  }

  shift=0.5*kshift;
  if(conservative) {
    ip_low=nint_he(nr*phi_low/(2*M_PI)-shift);
    ip_hi=nint_he(nr*phi_hi/(2*M_PI)-shift);
    ip_low=modu_he(ip_low,nr);
    ip_hi=modu_he(ip_hi,nr);
  }
  else {
    ip_low=(int)(nr*phi_low/(2*M_PI)-shift)+1;
    ip_hi=(int)(nr*phi_hi/(2*M_PI)-shift);
    diff=modu_he((ip_low-ip_hi),nr);
    if(diff<0) diff+=nr;
    if((diff==1)&&(dphi*nr<M_PI)) {
      *nir=0;
      return;
    }
    if(ip_low>=nr) ip_low-=nr;
    if(ip_hi<0) ip_hi+=nr;
  }

  if(ip_low>ip_hi) to_top=1;
  ip_low+=ipix1;
  ip_hi+=ipix1;

  if(to_top) {
    nir1=ipix2-ip_low+1;
    nir2=ip_hi-ipix1+1;
    (*nir)=nir1+nir2;

    if(*nir>nir_here)
      report_error(1,"Not enough memory in listir\n");
    for(jj=0;jj<nir1;jj++)
      listir[jj]=ip_low+jj;
    for(jj=nir1;jj<(*nir);jj++)
      listir[jj]=ipix1+jj-nir1;
  }
  else {
    (*nir)=ip_hi-ip_low+1;

    if(*nir>nir_here)
      report_error(1,"Not enough memory in listir\n");
    for(jj=0;jj<(*nir);jj++)
      listir[jj]=ip_low+jj;
  }

  return;
}

static double wrap_phi(double phi)
{
  if(phi>2*M_PI)
    return wrap_phi(phi-2*M_PI);
  else if(phi<0)
    return wrap_phi(phi+2*M_PI);
  else
    return phi;
}

static int he_ring_num(long nside,double z)
{
  //Returns ring index for normalized height z
  int iring;

  iring=(int)(nside*(2-1.5*z)+0.5);
  if(z>0.66666666) {
    iring=(int)(nside*sqrt(3*(1-z))+0.5);
    if(iring==0) iring=1;
  }

  if(z<-0.66666666) {
    iring=(int)(nside*sqrt(3*(1+z))+0.5);
    if(iring==0) iring=1;
    iring=4*nside-iring;
  }

  return iring;
}

static void he_query_disc(int nside,double cth0,double phi,double radius,
			  int *listtot,int *nlist,int inclusive)
{
  double phi0;
  int irmin,irmax,iz,nir,ilist;
  double radius_eff,a,b,c,cosang;
  double dth1,dth2,cosphi0,cosdphi,dphi;
  double rlat0,rlat1,rlat2,zmin,zmax,z;
  int *listir;

  phi0=wrap_phi(phi);
  listir=&(listtot[*nlist]);

  if((radius<0)||(radius>M_PI))
    report_error(1,"The angular radius is in RADIAN, and should lie in [0,M_PI]!");

  dth1=1/(3*((double)(nside*nside)));
  dth2=2/(3*((double)nside));

  if(inclusive)
    radius_eff=radius+1.071*M_PI/(4*((double)nside)); //TODO:check this
  else
    radius_eff=radius;
  cosang=cos(radius_eff);
  
  //Circle center
  cosphi0=cos(phi0);
  a=1-cth0*cth0;
  
  //Coord z of highest and lowest points in the disc
  rlat0=asin(cth0); //TODO:check
  rlat1=rlat0+radius_eff;
  rlat2=rlat0-radius_eff;

  if(rlat1>=0.5*M_PI)
    zmax=1;
  else
    zmax=sin(rlat1);
  irmin=he_ring_num(nside,zmax);
  if(irmin<2)
    irmin=1;
  else
    irmin--;

  if(rlat2<=-0.5*M_PI)
    zmin=-1;
  else
    zmin=sin(rlat2);
  irmax=he_ring_num(nside,zmin);
  if(irmax>(4*nside-2))
    irmax=4*nside-1;
  else
    irmax++;

  ilist=0;
  //Loop on ring number
  for(iz=irmin;iz<=irmax;iz++) {
    int kk;

    if(iz<=nside-1) //North polar cap
      z=1-iz*iz*dth1;
    else if(iz<=3*nside) //Tropical band + equator
      z=(2*nside-iz)*dth2;
    else
      z=-1+dth1*(4*nside-iz)*(4*nside-iz);
    
    //phi range in the disc for each z
    b=cosang-z*cth0;
    c=1-z*z;
    if(cth0==1) {
      dphi=M_PI;
      if(b>0) continue; //Out of the disc
    }
    else {
      cosdphi=b/sqrt(a*c);
      if(fabs(cosdphi)<=1)
	dphi=acos(cosdphi);
      else {
	if(cosphi0<cosdphi) continue; //Out of the disc
	dphi=M_PI; 
      }
    }

    //Find pixels in the disc
    nir=*nlist;
    he_in_ring(nside,iz,phi0,dphi,listir,&nir);
    //    printf("%lf %lf %d\n",dphi,cosdphi,nir);

    if(*nlist<ilist+nir) {
      report_error(1,"Not enough memory in listtot %d %d %lf %lf %lf %d\n",
		   *nlist,ilist+nir,radius,cth0,phi,nside);
    }
    for(kk=0;kk<nir;kk++) {
      listtot[ilist]=listir[kk];
      ilist++;
    }
  }

  *nlist=ilist;
}

static inline int th2bin_scaled(double cth,double scale,int logbin,
				double n_logint,double log_x_max,
				int nb_x,double i_x_max)
{
  int ix;
  cth=(MIN((1.),(cth)));
  
  if(logbin) {
    if(cth!=1) {
#ifdef _TRUE_ACOS
      cth=log10(acos((MIN(1,cth)))*scale);
#else //_TRUE_ACOS
      cth=1-MIN(1,cth);
      cth=0.5*log10((2*cth+0.3333333*cth*cth+
		     0.0888888889*cth*cth*cth)*scale*scale);
#endif //_TRUE_ACOS
      ix=(int)(n_logint*(cth-log_x_max)+nb_x);
    }
    else ix=-1;
  }
  else {
#ifdef _TRUE_ACOS
    cth=acos((MIN(1,cth)))*scale;
#else //_TRUE_ACOS
    cth=1-MIN(1,cth);
    cth=scale*sqrt(2*cth+0.333333333*cth*cth+
		   0.08888888889*cth*cth*cth);
#endif //_TRUE_ACOS
    ix=(int)(cth*nb_x*i_x_max);
  }
  
  return ix;
}

static inline double getrad(double cth,double scale)
{
  double r=(MIN((1.),(cth)));

#ifdef _TRUE_ACOS
  r=acos((MIN(1,r)))*scale;
#else //_TRUE_ACOS
  r=1-MIN(1,r);
  r=scale*sqrt(2*r+0.333333333*r*r+0.08888888889*r*r*r);
#endif //_TRUE_ACOS
  
  return r;
}

void cute_correlation_scaled(long ngal,double *pos,long nside,double *fld,double *msk,
			     double xmin,double xmax,int nx,int do_log,
			     double *hf_th,double *hm_th)
{
  long npix=he_nside2npix(nside);
  double *pos_pix=my_malloc(3*npix*sizeof(double));

#pragma omp parallel default(none) \
  shared(pos_pix,npix,nside)
  {
    long ip;

#pragma omp for
    for(ip=0;ip<npix;ip++) {
      double *v=&(pos_pix[3*ip]);
      he_pix2vec_ring(nside,ip,v);
    } //end omp for
  } //end omp parallel

  int ih;
  for(ih=0;ih<nx;ih++) {
    hf_th[ih]=0;
    hm_th[ih]=0;
  }

#pragma omp parallel default(none)			\
  shared(ngal,pos,fld,msk,nside,xmin,xmax,nx,do_log)	\
  shared(hf_th,hm_th,npix,pos_pix)
  {
    int i;
    double *hf_th_thr=my_calloc(nx,sizeof(double));
    double *hm_th_thr=my_calloc(nx,sizeof(double));
    int lenlist0=2*npix;
    int *listpix=my_malloc(lenlist0*sizeof(int));

    int logbin=do_log;
    double log_x_max=log10(xmax);
    double i_x_max=1./xmax;
    double n_logint=-1;
    if(do_log)
      n_logint=nx/log10(xmax/xmin);

#pragma omp for
    for(i=0;i<ngal;i++) {
      int j;
      double pos_g[3];
      int lenlist_half=lenlist0/2;
      double cth_g=pos[4*i+0];
      double phi_g=pos[4*i+1];
      double wei_g=pos[4*i+2];
      double sca_g=pos[4*i+3];
      double thmax_g=xmax/sca_g;
      double cthmax_g=cos(thmax_g);
      if(thmax_g*1.2>M_PI)
	printf("%lE %lE %lE %lE\n", sca_g, 1./sca_g, 180/(M_PI*sca_g), 180*thmax_g/M_PI);
      pos_g[0]=sqrt(1-cth_g*cth_g)*cos(phi_g);
      pos_g[1]=sqrt(1-cth_g*cth_g)*sin(phi_g);
      pos_g[2]=cth_g;
      he_query_disc(nside,cth_g,phi_g,1.2*thmax_g,listpix,&lenlist_half,1);
      //printf("HII %d\n", lenlist_half);
      for(j=0;j<lenlist_half;j++) {
	int ipx=listpix[j];
	double *pos_p=&(pos_pix[3*ipx]);
	double prod=pos_g[0]*pos_p[0]+pos_g[1]*pos_p[1]+pos_g[2]*pos_p[2];
	if(prod>cthmax_g) {
	  int ix=th2bin_scaled(prod,sca_g,logbin,n_logint,log_x_max,nx,i_x_max);
	  if((ix<nx)&&(ix>=0)) {
	    //printf("%d\n", ix);
	    hf_th_thr[ix]+=wei_g*fld[ipx];
	    hm_th_thr[ix]+=wei_g*msk[ipx];
	  }
	}
      }
      /*
      */
    } //end omp for

#pragma omp critical
    {
      for(i=0;i<nx;i++) {
	hf_th[i]+=hf_th_thr[i];
	hm_th[i]+=hm_th_thr[i];
      }
    } //end omp critical

    free(hf_th_thr);
    free(hm_th_thr);
    free(listpix);
  } //end omp parallel

  free(pos_pix);
}

void cute_correlation_scaled_2D(long ngal,double *pos,long nside,double *fld,double *msk,
				double xmin,double xmax,int nx,int na,int do_log,
				double *hf_th,double *hm_th)
{
  long npix=he_nside2npix(nside);
  double *pos_pix=my_malloc(4*npix*sizeof(double));

#pragma omp parallel default(none) \
  shared(pos_pix,npix,nside)
  {
    long ip;

#pragma omp for
    for(ip=0;ip<npix;ip++) {
      double *v=&(pos_pix[4*ip]);
      he_pix2vec_ring(nside,ip,v);
      v[3]=atan2(v[1],v[0]);
    } //end omp for
  } //end omp parallel

  int ih;
  for(ih=0;ih<nx*na;ih++) {
    hf_th[ih]=0;
    hm_th[ih]=0;
  }

#pragma omp parallel default(none)			\
  shared(ngal,pos,fld,msk,nside,xmin,xmax,nx,do_log)	\
  shared(hf_th,hm_th,npix,pos_pix,na)
  {
    int i;
    double *hf_th_thr=my_calloc(nx*na,sizeof(double));
    double *hm_th_thr=my_calloc(nx*na,sizeof(double));
    int lenlist0=npix/2;
    int *listpix=my_malloc(lenlist0*sizeof(int));

    int logbin=do_log;
    double log_x_max=log10(xmax);
    double i_x_max=1./xmax;
    double n_logint=-1;
    if(do_log)
      n_logint=nx/log10(xmax/xmin);

#pragma omp for
    for(i=0;i<ngal;i++) {
      int j;
      double pos_g[3];
      int lenlist_half=lenlist0/2;
      double cth_g=pos[4*i+0];
      double sth_g=sqrt(1-cth_g*cth_g);
      double phi_g=pos[4*i+1];
      double wei_g=pos[4*i+2];
      double sca_g=pos[4*i+3];
      double thmax_g=xmax/sca_g;
      double cthmax_g=cos(thmax_g);
      pos_g[0]=sth_g*cos(phi_g);
      pos_g[1]=sth_g*sin(phi_g);
      pos_g[2]=cth_g;
      he_query_disc(nside,cth_g,phi_g,1.2*thmax_g,listpix,&lenlist_half,1);
      for(j=0;j<lenlist_half;j++) {
	int ipx=listpix[j];
	double *pos_p=&(pos_pix[4*ipx]);
	double cth_gf=MIN(1,(pos_g[0]*pos_p[0]+pos_g[1]*pos_p[1]+pos_g[2]*pos_p[2]));
	double cth_f=pos_p[2];
	if(cth_gf>cthmax_g) {
	  int ix=th2bin_scaled(cth_gf,sca_g,logbin,n_logint,log_x_max,nx,i_x_max);	
	  if((ix<nx)&&(ix>=0)) {
	    int ia;
	    double alpha=0;
	    if((sth_g==0) || (cth_gf>=1) || (cth_gf<=-1))
	      ia=0;
	    else {
	      double sth_a;
	      double phi_f=pos_p[3];
	      double sth_gf=sqrt(1-cth_gf*cth_gf);
	      double cth_a=(cth_f-cth_gf*cth_g)/(sth_g*sth_gf);
	      if((cth_a>=1) || (cth_a<=-1))
		sth_a=0;
	      else {
		if(sin(phi_g-phi_f)>0)
		  sth_a= sqrt(1-cth_a*cth_a);
		else
		  sth_a=-sqrt(1-cth_a*cth_a);
	      }
	      alpha=atan2(sth_a,cth_a);
	      if(alpha<0)
		alpha+=2*M_PI;
	      else if(alpha>=2*M_PI)
		alpha-=2*M_PI;
	      ia=(int)(na*alpha/(2*M_PI));
	      if(ia>=na)
		ia-=na;
	    }
	    hf_th_thr[na*ix+ia]+=wei_g*fld[ipx];
	    hm_th_thr[na*ix+ia]+=wei_g*msk[ipx];
	  }
	}
      }
    } //end omp for

#pragma omp critical
    {
      for(i=0;i<nx*na;i++) {
	hf_th[i]+=hf_th_thr[i];
	hm_th[i]+=hm_th_thr[i];
      }
    } //end omp critical

    free(hf_th_thr);
    free(hm_th_thr);
    free(listpix);
  } //end omp parallel

  free(pos_pix);
}

void cute_line_correlation(long nside, double **fld, double *msk,
			   double xmin, double xmax, int nx, int do_log,
			   double *hf_th, double *hm_th, int per_bin)
{
  long npix=he_nside2npix(nside);
  double *pos_pix=my_malloc(3*npix*sizeof(double));
  double *cth_pix=my_malloc(npix*sizeof(double));
  double *phi_pix=my_malloc(npix*sizeof(double));

  //Calculate coordinates of all pixels in map
#pragma omp parallel default(none)		\
  shared(pos_pix,cth_pix,phi_pix,npix,nside)
  {
    long ip;

#pragma omp for
    for(ip=0;ip<npix;ip++) {
      double *v=&(pos_pix[3*ip]);
      he_pix2vec_ring(nside,ip,v);
      cth_pix[ip]=v[2];
      phi_pix[ip]=atan2(v[1],v[0]);
    } //end omp for
  } //end omp parallel

  //Zero all histograms
  int ih;
  for(ih=0;ih<nx;ih++) {
    hf_th[ih]=0;
    hm_th[ih]=0;
  }

  //Parallel region
#pragma omp parallel default(none)			\
  shared(fld,msk,nside,npix,xmin,xmax,nx,do_log)	\
  shared(hf_th,hm_th,pos_pix,cth_pix,phi_pix,per_bin)
  {
    int i1;
    double *hf_th_thr=my_calloc(nx,sizeof(double));
    double *hm_th_thr=my_calloc(nx,sizeof(double));
    int lenlist0=2*npix;
    int *listpix=my_malloc(lenlist0*sizeof(int));

    int logbin=do_log;
    double log_x_max=log10(xmax);
    double i_x_max=1./xmax;
    double n_logint=-1;
    if(do_log)
      n_logint=nx/log10(xmax/xmin);

#pragma omp for schedule(dynamic)
    for(i1=0;i1<npix;i1++) {
      if(msk[i1]<=0)
	continue;
      int j, lenlist_half=lenlist0/2;
      double *pos1=&(pos_pix[3*i1]);
      he_query_disc(nside,cth_pix[i1],phi_pix[i1],2.2*xmax,listpix,&lenlist_half,1);
      for(j=0;j<lenlist_half;j++) {
	long i2;
	double *f;
	double pos2[3];
	int ax,ix,i3=listpix[j];
	double *pos3=&(pos_pix[3*i3]);
	double prod=0;
	if(i3<=i1)
	  continue;
	if(msk[i3]<=0)
	  continue;
	for(ax=0;ax<3;ax++) {
	  pos2[ax]=(pos1[ax]+pos3[ax])*0.5;
	  prod+=pos1[ax]*pos3[ax];
	}
	vec2pix_ring(nside,pos2,&i2);
	if(msk[i2]<=0)
	  continue;
	ix=th2bin_scaled(prod,0.5,logbin,n_logint,log_x_max,nx,i_x_max);
	if((ix<nx) && (ix>=0)) {
	  if(per_bin)
	    f=fld[ix];
	  else
	    f=fld[0];
	  hf_th_thr[ix]+=f[i1]*f[i2]*f[i3];
	  hm_th_thr[ix]+=msk[i1]*msk[i2]*msk[i3];
	}
      }
    } //end omp for

#pragma omp critical
    {
      for(i1=0;i1<nx;i1++) {
	hf_th[i1]+=hf_th_thr[i1];
	hm_th[i1]+=hm_th_thr[i1];
      }
    } //end omp critical

    free(hf_th_thr);
    free(hm_th_thr);
    free(listpix);
  } //end omp parallel

  free(pos_pix);
  free(cth_pix);
  free(phi_pix);
}

void cute_correlation(long nside, double *fld, double *msk,
		      double xmin, double xmax, int nx, int do_log,
		      double *hf_th, double *hm_th)
{
  long npix=he_nside2npix(nside);
  double *pos_pix=my_malloc(3*npix*sizeof(double));
  double *cth_pix=my_malloc(npix*sizeof(double));
  double *phi_pix=my_malloc(npix*sizeof(double));

  //Calculate coordinates of all pixels in map
#pragma omp parallel default(none)		\
  shared(pos_pix,cth_pix,phi_pix,npix,nside)
  {
    long ip;

#pragma omp for
    for(ip=0;ip<npix;ip++) {
      double *v=&(pos_pix[3*ip]);
      he_pix2vec_ring(nside,ip,v);
      cth_pix[ip]=v[2];
      phi_pix[ip]=atan2(v[1],v[0]);
    } //end omp for
  } //end omp parallel

  //Zero all histograms
  int ih;
  for(ih=0;ih<nx;ih++) {
    hf_th[ih]=0;
    hm_th[ih]=0;
  }

  //Parallel region
#pragma omp parallel default(none)			\
  shared(fld,msk,nside,npix,xmin,xmax,nx,do_log)	\
  shared(hf_th,hm_th,pos_pix,cth_pix,phi_pix)
  {
    int i1;
    double *hf_th_thr=my_calloc(nx,sizeof(double));
    double *hm_th_thr=my_calloc(nx,sizeof(double));
    int lenlist0=2*npix;
    int *listpix=my_malloc(lenlist0*sizeof(int));

    int logbin=do_log;
    double log_x_max=log10(xmax);
    double i_x_max=1./xmax;
    double n_logint=-1;
    if(do_log)
      n_logint=nx/log10(xmax/xmin);

#pragma omp for schedule(dynamic)
    for(i1=0;i1<npix;i1++) {
      if(msk[i1]<=0)
	continue;
      int j, lenlist_half=lenlist0/2;
      double *pos1=&(pos_pix[3*i1]);
      he_query_disc(nside,cth_pix[i1],phi_pix[i1],1.2*xmax,listpix,&lenlist_half,1);
      for(j=0;j<lenlist_half;j++) {
	int ax,ix,i2=listpix[j];
	double *pos2=&(pos_pix[3*i2]);
	double prod=0;
	if(i2<=i1)
	  continue;
	if(msk[i2]<=0)
	  continue;
	for(ax=0;ax<3;ax++)
	  prod+=pos1[ax]*pos2[ax];
	ix=th2bin_scaled(prod,1.0,logbin,n_logint,log_x_max,nx,i_x_max);
	if((ix<nx) && (ix>=0)) {
	  hf_th_thr[ix]+=fld[i1]*fld[i2];
	  hm_th_thr[ix]+=msk[i1]*msk[i2];
	}
      }
    } //end omp for

#pragma omp critical
    {
      for(i1=0;i1<nx;i1++) {
	hf_th[i1]+=hf_th_thr[i1];
	hm_th[i1]+=hm_th_thr[i1];
      }
    } //end omp critical

    free(hf_th_thr);
    free(hm_th_thr);
    free(listpix);
  } //end omp parallel

  free(pos_pix);
  free(cth_pix);
  free(phi_pix);
}

/*
static void get_ring_limits(long nside,int iz,long *ip_lo,long *ip_hi)
{
  long ir;
  long ipix1,ipix2;
  long npix=12*nside*nside;
  long ncap=2*nside*(nside-1);

  if((iz>=nside)&&(iz<=3*nside)) { //eqt
    ir=iz-nside+1;
    ipix1=ncap+4*nside*(ir-1);
    ipix2=ipix1+4*nside-1;
  }
  else {
    if(iz<nside) { //north
      ir=iz;
      ipix1=2*ir*(ir-1);
      ipix2=ipix1+4*ir-1;
    }
    else { //south
      ir=4*nside-iz;
      ipix1=npix-2*ir*(ir+1);
      ipix2=ipix1+4*ir-1;
    }
  }

  *ip_lo=ipix1;
  *ip_hi=ipix2;
}
*/

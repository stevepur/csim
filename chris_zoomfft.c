
Since the chirp-z transform (czt) gets called all the time, we set up work space and
fftw plans once, up front, based on the size of the problem:

void czt_setup(int m, int k)
{
 int n,nfft;

 if (m>=k) n=m; else n=k;
 nfft = pow(2,ceil(log2(m+k-1)));  // next largest power of 2, for ffts
 #pragma omp parallel
 {
#pragma omp critical   // let's be safe
   {
     fy = malloc(nfft*sizeof(complex));
     fv = malloc(nfft*sizeof(complex));
     yy = malloc(nfft*sizeof(complex));
     iww = malloc(nfft*sizeof(complex));
         ww = malloc((m+n-1)*sizeof(complex));
     g = malloc((m+n-1)*sizeof(complex));
         aa = malloc(m*sizeof(complex));
   }
 }

 // single shared copy of plans, for use with new-array execute
 p1 = fftw_plan_dft_1d(nfft,yy,fy,FFTW_FORWARD,FFTW_ESTIMATE);
 p2 = fftw_plan_dft_1d(nfft,iww,fv,FFTW_FORWARD,FFTW_ESTIMATE);
 p3 = fftw_plan_dft_1d(nfft,fy,g,FFTW_BACKWARD,FFTW_ESTIMATE);

}



The czt is mainly used for a "zoomFFT", which is how it usually gets called at
the top level. N is array size, set elsewhere - typically 1024.

void zoomFFT_realunits(double *x, double *y, void *in, double *xi, double *eta, void *out,
                      double f, double lambda)
{
 int i,j;
 double dxi = xi[1]-xi[0];
 double deta = eta[1]-eta[0];
 double dx = x[1]-x[0];
 double dy = y[1]-y[0];

 complex Ax = cexp(2*pi*I * dx *xi[0]/(lambda*f));
 complex Wx = cexp(-2*pi*I * dx *dxi/(lambda*f));
 complex Ay = cexp(2*pi*I * dy *eta[0]/(lambda*f));
 complex Wy = cexp(-2*pi*I * dy *deta/(lambda*f));

 complex (*outarr)[N] = (complex(*)[N])out;

 czt2d(in,N,N,Wx,Ax,Wy,Ay,out);
#pragma omp parallel for private(j) collapse(2)
 for (i=0;i<N;++i)
   for (j=0;j<N;++j)
     outarr[i][j] *= cexp(-2*pi*I * (x[0]*xi[i]+y[0]*eta[j])/(lambda*f))
       * dx*dy/(I*lambda*f);
}



Like a 2D FFT, czt2d does 1D czts in memory direction, then transposes for a second
round of 1D czts:

int czt2d(complex *x, int m, int k, complex wx, complex ax, complex wy, complex ay, complex *out)
{
 // symmetric case: square input, same {k,w,a} in both dims
 // trivial mod for different {m,k,w,a} in each dim
 int i,j;
 complex (*in)[m] = (complex(*)[m])x;
 complex (*o)[m] = (complex(*)[m])out;
 complex (*tmp)[m] = (complex(*)[m])malloc(m*m*sizeof(complex));

#pragma omp parallel
 {
   complex tv[m];
#pragma omp for
   for (i=0;i<m;++i)
     czt(in[i],m,k,wy,ay,tmp[i]);
   #pragma omp for private (j)
   for (i=0;i<m;++i){
     for (j=0;j<m;++j)
       tv[j] = tmp[j][i];
     czt(tv,m,k,wx,ax,o[i]);
   }
 }

 free(tmp);

 return 0;
}


The real business is in the 1D czt:

void czt(complex *x, int m, int k, complex w, complex a, complex *out)
{
 int i,j,n,nfft;
 double kk2;

 if (m>=k) n=m; else n=k;
 nfft = pow(2,ceil(log2(m+k-1)));  // next largest power of 2, for ffts

 /////

 for (i=0;i<m+n-1;++i){
   kk2 = pow(i-m+1,2)/2.;
   ww[i] = cpow(w,kk2);
 }

 for (i=0;i<m;++i)
   aa[i] = cpow(a,-i)*ww[m-1+i];

 memset(yy,0,nfft*sizeof(complex));  // zero pad
   for (i=0;i<m;++i)
   yy[i] = x[i]*aa[i];

 fftw_execute_dft(p1,yy,fy);

 memset(iww,0,nfft*sizeof(complex));  // zero pad

 for (i=0;i<m+n-1;++i)
   iww[i] = 1./ww[i];

 fftw_execute_dft(p2,iww,fv);

 for (i=0;i<nfft;++i)
   fy[i]*=fv[i];

 fftw_execute_dft(p3,fy,g);

 for (i=0;i<k;++i)
   g[m-1+i]*=ww[m-1+i]/nfft;

 for (i=0;i<k;++i)
   out[i] = g[m-1+i];

}


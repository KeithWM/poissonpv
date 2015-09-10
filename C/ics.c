#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mat.h"
#include "normal.h"

#define TARGETM 0
#define TARGETN 0
#define TARGETMN 1
#define GAMMA 1

void writeArrayToMatlab(MATFile *pmatfile, double *store, char *name, int m, int n)
{
  mxArray *matlabArray;
  
  matlabArray = (mxArray *) mxCreateDoubleMatrix(m,n,mxREAL); /* creates a mxn MATLAB array of reals */
  memcpy((void *)(mxGetPr(matlabArray)),   (void *)store, m*n*sizeof(double)); /* copies the memory (m*n) at Q to tseries */
  matPutVariable(pmatfile,name,matlabArray); /* puts the array matlabx into the .mat file as "dt" */
  mxDestroyArray(matlabArray); /* clears the array (for memory I presume) */  
}

void conserved(double *x, double *y, double *z, double *Gamma, double *H, double *J, int N)
{
  /* computes the conserved quantities H, J1 and J2 */
  /* IGNORES OMEGA! */
  double pi = acos(-1);
  int i, j;
  
  H[0] = 0.;
  J[0] = 0.;
  J[1] = 0.;
  J[2] = 0.;
  
  for(i=0; i<N; ++i)
  {
    for(j=0; j<i; ++j) {
      H[0]   -= Gamma[i]*Gamma[j]*log( 2 * ( 1. - x[i]*x[j] - y[i]*y[j] - z[i]*z[j] ) );
    }
    J[0] += Gamma[i]*x[i];
    J[1] += Gamma[i]*y[i];
    J[2] += Gamma[i]*z[i];
  }
  
  H[0] /= 4*pi;
}

int Jcloser(double *J, double *Jpartold, double *Jpart, double Jtarget)
{
  double distnew = 0., distold = 0., tmp;
  int i;
  
  for(i=0;i<3;++i)
  {
    if(i==2)
      tmp = J[i] - Jtarget;
    else
      tmp = J[i];
    distold += tmp*tmp;
    tmp+= -Jpartold[i] + Jpart[i];
    distnew += tmp*tmp;
  }  
  
  return (distnew < distold);
}
//       if(fabs(H[0] - Hpartold + Hpart - Etarget) < fabs(H[0]-Etarget)) {

void Hpartfun(double *x, double *y, double *z, double *Gamma, int k, int l, double *H, int N)
{
  /* computes the part of H due to pair kl */
  double pi = acos(-1);
  int i;
  
  H[0] = 0.;
  
  for(i=0; i<N; ++i)
  {
    if((i != k) && (i != l)) {
      H[0]   -= Gamma[i]*Gamma[k]*log( 2 * ( 1. - x[i]*x[k] - y[i]*y[k] - z[i]*z[k] ) );
      H[0]   -= Gamma[i]*Gamma[l]*log( 2 * ( 1. - x[i]*x[l] - y[i]*y[l] - z[i]*z[l] ) );
    }
  }
  
  H[0] /= 4*pi;
}

void Jpartfun(double *x, double *y, double *z, double *Gamma, int k, int l, double *J)
{  
  J[0] = x[k]*Gamma[k] + x[l]*Gamma[l];
  J[1] = y[k]*Gamma[k] + y[l]*Gamma[l];
  J[2] = z[k]*Gamma[k] + z[l]*Gamma[l];
}

void Hrotate(double *x, double *y, double *z, double *Gamma, int i, int j, double theta)
// rotate two point vortices while maintaining J_{ij}
{
  double cx, cy, cz, cnorm, s, c, cc1, cc2, xtemp, ytemp;
  
  cx = theta * (Gamma[i]*x[i] + Gamma[j]*x[j]);
  cy = theta * (Gamma[i]*y[i] + Gamma[j]*y[j]);
  cz = theta * (Gamma[i]*z[i] + Gamma[j]*z[j]);
  
//printf("c[%d]: (%f+%f,%f+%f,%f+%f) = (%f,%f,%f)\n", i, Gamma[i]*x[i], Gamma[j]*x[j], Gamma[i]*y[i], Gamma[j]*y[j], Gamma[i]*z[i], Gamma[j]*z[j], cx, cy, cz);
  
  cnorm = sqrt(cx*cx + cy*cy + cz*cz);
  if(cnorm != 0.) {
    cx = cx/cnorm;
    cy = cy/cnorm;
    cz = cz/cnorm;
    
    s = sin(cnorm);
    c = 1.-cos(cnorm);
    
    cc1 = cx*x[i] + cy*y[i] + cz*z[i];
    cc2 = cx*x[j] + cy*y[j] + cz*z[j];
    
    xtemp = x[i] + s * ( cy*z[i]-cz*y[i] ) - c * ( x[i]-cc1*cx );
    ytemp = y[i] + s * ( cz*x[i]-cx*z[i] ) - c * ( y[i]-cc1*cy );
    z[i]  = z[i] + s * ( cx*y[i]-cy*x[i] ) - c * ( z[i]-cc1*cz );
    x[i]  = xtemp;
    y[i]  = ytemp;
    
    xtemp = x[j] + s * ( cy*z[j]-cz*y[j] ) - c * ( x[j]-cc2*cx );
    ytemp = y[j] + s * ( cz*x[j]-cx*z[j] ) - c * ( y[j]-cc2*cy );
    z[j]  = z[j] + s * ( cx*y[j]-cy*x[j] ) - c * ( z[j]-cc2*cz );
    x[j]  = xtemp;
    y[j]  = ytemp;
  } else {
    printf("cnorm zero for pair (%d,%d)\n", i, j);
  }
}

double mindistance(double *x, double *y, double *z, int M, int N) {
  double d, min = 10.;
  int i,j;
  
  for (i=M; i<N; ++i) {
    for (j=M; j<i; ++j) {
      d = 1 - x[i]*x[j]-y[i]*y[j]-z[i]*z[j];
      min = (min < d) ? min : d;
    }
  }
  return min;
}

void HiterPlace(double *x, double *y, double *z, double *Gamma, double *H, double *J, double Etarget, int M, int N, int *seed)
// rotate pairs such that the Hamiltonian appraoches the target
{
  int h, i, j, j0, k=0, hmax = 10;
  double *theta, *is, pi=acos(-1), r=1;
  double u1, v1, w1, u2, v2, w2, Hpart, Hpartold;
  
  theta = r8vec_uniform_01_new(hmax*hmax, seed);
  is = r8vec_uniform_01_new(N-M, seed);
  
  for (h=0; h<hmax; ++h) {
//    printf("%d out of %d, i: %d, H: %e\n", h, N-M, i, H[0]);
    i = round(is[h]*(N-M)-.5) + M;
    j0 = (M > i-hmax) ? M : i-hmax;
    for (j=j0; j<i; ++j) {
      theta[k] = 2*pi*theta[k];
      u1 = x[i]; v1 = y[i]; w1 = z[i];
      u2 = x[j]; v2 = y[j]; w2 = z[j];
      
      Hpartfun(x, y, z, Gamma, i, j, &Hpartold, N);
      Hrotate(x, y, z, Gamma, i, j, theta[k]);
      Hpartfun(x, y, z, Gamma, i, j, &Hpart, N);
      
      if(fabs(H[0] - Hpartold + Hpart - Etarget) < fabs(H[0]-Etarget)) {
        /* accept rotation */
        //printf("ACCEPTING (%d,%d)\n", i,j);
        H[0] += -Hpartold + Hpart;
        r = H[0] - Etarget;
        //printf("r: %.2e\n", r);
      } else {
        /* reject rotation */
        //printf("rejecting (%d,%d)\n", i,j);
        x[i] = u1; y[i] = v1; z[i] = w1;
        x[j] = u2; y[j] = v2; z[j] = w2;
      }
      ++k;
    }
  }
  free(theta);
  free(is);
  printf("Round end! (H: %3.3e, H-Etarget: %3.3e) for %d. Mindistance: %e\n", H[0], H[0]-Etarget, N, mindistance(x,y,z,M,N));
}

void HPlace(double *x, double *y, double *z, double *Gamma, double *H, double *J, double Etarget, double tol, int M, int N, int *seed)
// rotate pairs such that the Hamiltonian appraoches the target
{
  int i, j, k, l=0;
  double *theta, pi=acos(-1), r=1;
  double *u, *v, *w, Hold;
  
  u     = (double *) malloc(N*(sizeof(double)));
  v     = (double *) malloc(N*(sizeof(double)));
  w     = (double *) malloc(N*(sizeof(double)));
  
  conserved(x, y, z, Gamma, H, J, N);
  Hold = *H;
  while(fabs(*H-Etarget) > tol && l < 3000) {
    k = 0;
    theta = r8vec_uniform_01_new((N-M)*(N-M-1)/2, seed);
    
    for (i=M; i<N; ++i) {
      u[i] = x[i];
      v[i] = y[i];
      w[i] = z[i];
//printf("x[%d]: (%f,%f,%f)\n", i, x[i], y[i], z[i]);
    }
    for (i=M; i<N; ++i) {
      for (j=M; j<i; ++j) {
        theta[k] = 2*pi*r*theta[k];
//printf("r: %f, theta[%d]: %f\n", r, k, theta[k]);
        Hrotate(x, y, z, Gamma, i, j, theta[k]);
        ++k;
      }
    }
    conserved(x, y, z, Gamma, H, J, N);
    //printf("Round test %d! (H: %3.3e, H-Etarget: %3.3e) for %d. Mindistance: %f\n", l, H[0], H[0]-Etarget, N, mindistance(x,y,z,M,N));
    if(fabs(*H - Etarget) < fabs(Hold - Etarget)) {
      /* accept rotation */
      //printf("ACCEPTING\n");
      Hold = *H;
      r = .01*(Hold - Etarget);
      //printf("r: %f\n", r);
    } else {
      /* reject rotation */
      //printf("rejecting\n");
      for (i=M; i<N; ++i) {
        x[i] = u[i];
        y[i] = v[i];
        z[i] = w[i];
//printf("x[%d]: (%f,%f,%f)\n", i, x[i], y[i], z[i]);
      }
    }
    conserved(x, y, z, Gamma, H, J, N);
    printf("Round end %d! (H: %3.3e, H-Etarget: %3.3e) for %d. Mindistance: %f\n", l++, H[0], H[0]-Etarget, N, mindistance(x,y,z,M,N));
  }
}

void JiterPlace(double *x, double *y, double *z, double *Gamma, double *H, double *J, double Jtarget, int M, int N, int *seed)
// rotate pairs such that the Hamiltonian appraoches the target
{
  int i, j, k=0;
  double pi=acos(-1);
  double u1, v1, w1, unorm, Jres[3], Jt[3];
  
  conserved(x, y, z, Gamma, H, J, N);
  Jres[0] = -J[0];
  Jres[1] = -J[1];
  Jres[2] = Jtarget-J[2];
  
  for(i=M; i<N; ++i) {
    u1 = Gamma[i]*x[i] + Jres[0];
    v1 = Gamma[i]*y[i] + Jres[1];
    w1 = Gamma[i]*z[i] + Jres[2];
    u1/= Gamma[i];
    v1/= Gamma[i];
    w1/= Gamma[i];
    
    unorm = sqrt(u1*u1 + v1*v1 + w1*w1);
//printf("x[%d]: (%f,%f,%f)\n", i, x[i], y[i], z[i]);
// printf("Jres: (%3.3e,%3.3e,%3.3e),\n", Jres[0],Jres[1],Jres[2]);
// printf("u[%d]: (%f,%f,%f), |u[%d]|: %f\n", i, u1, v1, w1, unorm);
    
    u1/= unorm;
    v1/= unorm;
    w1/= unorm;
    
    Jt[0] = Jres[0] + Gamma[i]*(x[i]-u1);
    Jt[1] = Jres[1] + Gamma[i]*(y[i]-v1);
    Jt[2] = Jres[2] + Gamma[i]*(z[i]-w1);
    
    if(Jt[0]*Jt[0] + Jt[1]*Jt[1] + Jt[2]*Jt[2] < Jres[0]*Jres[0] + Jres[1]*Jres[1] + Jres[2]*Jres[2]) {
      x[i] = u1;
      y[i] = v1;
      z[i] = w1;
      
      Jres[0] = Jt[0];
      Jres[1] = Jt[1];
      Jres[2] = Jt[2];
//printf("Jres: (%3.3e,%3.3e,%3.3e),\n", Jres[0],Jres[1],Jres[2]);
    }
//printf("x[%d]: (%f,%f,%f)\n", i, x[i], y[i], z[i]);
// printf("x[%d]: (%f,%f,%f)\n", i, x[i], y[i], z[i]);
  }
  J[0] = -Jres[0];
  J[1] = -Jres[1];
  J[2] = Jtarget-Jres[2];
  printf("Round end! (J: (%3.3e,%3.3e,%3.3e), J[2]-Jtarget: %3.3e) for %d\n", J[0],J[1],J[2], J[2]-Jtarget, N);
}

void perturbPlace(double *x, double *y, double *z, double *Gamma, double r, double *H, int M, int N, int *seed)
// perturb the initial placement, as it has all vortices directly opposite one another
{
  int i, j, jmin, k=0, jp = 10;
  double *theta, pi=acos(-1);
  
  for (i=M; i<N; ++i) {
    k = 0;
    theta = r8vec_uniform_01_new(jp, seed);
    jmin = (i > jp + M) ? i-jp : M;
    for (j=jmin; j<i; ++j) {
//      printf("Perturbing (%d,%d)\n", i, j);
      theta[k] = 2*pi*r*theta[k];      
//       printf("x[%d]: (%f,%f,%f) |x|-1: %f\t", i, x[i], y[i], z[i], x[i]*x[i] + y[i]*y[i] + z[i]*z[i] - 1);
//       printf("x[%d]: (%f,%f,%f) |x|-1: %f\n", j, x[j], y[j], z[j], x[j]*x[j] + y[j]*y[j] + z[j]*z[j] - 1);
//       printf("Rotating pair (%d,%d) by %f\n", i,j, theta[k]);
      Hrotate(x, y, z, Gamma, i, j, theta[k]);
//       printf("x[%d]: (%f,%f,%f) |x|-1: %f\t", i, x[i], y[i], z[i], x[i]*x[i] + y[i]*y[i] + z[i]*z[i] - 1);
//       printf("x[%d]: (%f,%f,%f) |x|-1: %f\n", j, x[j], y[j], z[j], x[j]*x[j] + y[j]*y[j] + z[j]*z[j] - 1);
      ++k;
    }
  }
  free(theta);
}

void initialPlace(double *x, double *y, double *z, int M, int N, double Jtarget, double *Gamma, int *pidum)
{
  int i, j=0, k;
  double *theta, pi=acos(-1), zmin, zmax;
  
#if GAMMA == 1
  theta = r8vec_uniform_01_new(2*(N-M), pidum); // a random vector, used later as BOTH theta AND z!
  
  if(Jtarget == 0.) { // place all vortices anipedian
    printf("Placing all vortices antipedian\n");
    for (i=0; i<(N-M)/2; ++i) {
      theta[j]   = 2*pi*theta[j];  // theta "proper"
      theta[j+1] = 1-2*theta[j+1]; // z
      x[M+2*i]   = sin(theta[j])*sqrt(1-theta[j+1]*theta[j+1]);
      y[M+2*i]   = cos(theta[j])*sqrt(1-theta[j+1]*theta[j+1]);
      z[M+2*i]   =               theta[j+1];
      x[M+2*i+1] = -x[M+2*i];
      y[M+2*i+1] = -y[M+2*i];
      z[M+2*i+1] = -z[M+2*i];
      j+= 2;
      //printf("x[%d]: (%f,%f,%f), x[%d]: (%f,%f,%f)\n", M+2*i, x[M+2*i], y[M+2*i], z[M+2*i], M+2*i+1, x[M+2*i+1], y[M+2*i+1], z[M+2*i+1]);
//       k = M+2*i;
//       printf("x[%d]: (%f,%f,%f) |x|-1: %f\n", k, x[k], y[k], z[k], x[k]*x[k] + y[k]*y[k] + z[k]*z[k] - 1);
//       k = M+2*i+1;
//       printf("x[%d]: (%f,%f,%f) |x|-1: %f\n", k, x[k], y[k], z[k], x[k]*x[k] + y[k]*y[k] + z[k]*z[k] - 1);
    } 
  } else { //place all vortices with J_{i,i+1} = 2*Jtarget/N, i.e. x_{i+1} = (2*Jt/N - G_i x_i)/G_{i+1}
    printf("Placing all vortices with partwise momentum\n");
    Jtarget*= 2./(N-M);
    for (i=0; i<(N-M)/2; ++i) {
      theta[i]   = 2*pi*theta[i];  
      z[M+2*i]   = Gamma[M+2*i]*Gamma[M+2*i+1] + Jtarget*Jtarget - Gamma[M+2*i+1]*Gamma[M+2*i+1];
      z[M+2*i]  /= 2*Jtarget*Gamma[M+2*i];
      x[M+2*i]   = sin(theta[i])*sqrt(1-z[M+2*i]*z[M+2*i]);
      y[M+2*i]   = cos(theta[i])*sqrt(1-z[M+2*i]*z[M+2*i]);
      x[M+2*i+1] = (0       - Gamma[M+2*i]*x[M+2*i])/Gamma[M+2*i+1];
      y[M+2*i+1] = (0       - Gamma[M+2*i]*y[M+2*i])/Gamma[M+2*i+1];
      z[M+2*i+1] = (Jtarget - Gamma[M+2*i]*z[M+2*i])/Gamma[M+2*i+1];
    }    
  }
#elif GAMMA > 1
//printf("|x[%d]: %f, |x[%d]|: %f\n", M+2*i, x[M+2*i]*x[M+2*i] + y[M+2*i]*y[M+2*i] + z[M+2*i]*z[M+2*i], M+2*i+1, x[M+2*i+1]*x[M+2*i+1] + y[M+2*i+1]*y[M+2*i+1] + z[M+2*i+1]*z[M+2*i+1]);
// printf("x[%d]: (%f,%f,%f), x[%d]: (%f,%f,%f)\n", M+2*i, x[M+2*i], y[M+2*i], z[M+2*i], M+2*i+1, x[M+2*i+1], y[M+2*i+1], z[M+2*i+1]);
//printf("Dx: (%f,%f,%f)\n", x[M+2*i] - x[M+2*i+1], y[M+2*i+1] - y[M+2*i], z[M+2*i] - z[M+2*i+1]);
  theta = r8vec_uniform_01_new(2*(N-M), pidum); // a random vector, used later as BOTH theta AND z!
    for (i=0; i<(N-M); ++i) {
      theta[2*i]   = 2*pi*theta[2*i];  // theta "proper"
      theta[2*i+1] = 1- 2*theta[2*i+1]; // z
      z[i]   = theta[2*i+1];
      x[i]   = sin(theta[2*i])*sqrt(1-z[i]*z[i]);
      y[i]   = cos(theta[2*i])*sqrt(1-z[i]*z[i]);
    }
#endif
}

void skewGamma(double *Gamma, double *x, double *y, double lambda, int N, int *seed)
{
  int i,j;
  double b,d, mean,std, pi = acos(-1);
    
  if(lambda == 0.) {
    r8vec_normal_01(N, seed, Gamma);
  } else {
    d = lambda/sqrt(1+lambda*lambda);
    b = sqrt(2/pi);

    mean = b*d;
    std  = sqrt(1-mean*mean);
    
    printf("Generating skew circulations, with lambda %f. Expected mean %f substracted and divided by expected variance %f\n", lambda, mean, std);
    
    i = j = 0;
    while(i<N) {
      if(j==0) {
        r8vec_normal_01(N, seed, x); // temporarily overload x
        r8vec_normal_01(N, seed, y); // temporarily overload y
      }
      x[j]/= lambda;
      if(x[j]<=lambda*y[j])
        Gamma[i++] = (y[j]-mean)/std;
      j=(j+1)%N;
    }  
  }
}

int readTargets(double *pE, double *pJ, char *energy)
{
  char path[50];
  const char *name;
  mxArray *pa;
  double *c;
  MATFile *pmat;
  
  sprintf(path, "/ufs/keith/Code/sphere/matlabdata/Lagrange_%s_M_MC100000_K5.mat", energy);
  printf("Reading file %s...\n\n", path);
  pmat = matOpen(path, "r");
  if (pmat == NULL) {
    printf("Error reopening file %s\n", path);
    return(1);
  }

  /* Get headers of all variables */
  printf("\nExamining the header for each variable:\n");
  while( (pa = matGetNextVariable(pmat, &name)) != NULL )
  {
//     printf("%s\n", name);
    if( strcmp(name, "c") == 0) {
      c = mxGetPr(pa);
    }
  }
  //printf("c: %f,%f,%f,%f,%f\n", c[0], c[1], c[2], c[3], c[4]);
  
  *pE = c[0];
  *pJ = sqrt(c[1]);

  if (matClose(pmat) != 0) {
    printf("Error closing file %s\n",path);
    return(1);
  }
  printf("Done reading array\n");
  return(0);
}

int main(int argc, char *argv[])
{
    /* constants */
  double pi = acos(-1);
  /* parameters */
  int M = 8, N = 8, L0=0, L = 1;
  double EMtarget = 0., ENtarget = 0., Jtarget = 0., tol = 1.e-10, lambda = 2.;
  double GammaA = 1.;
  double GammaB = 1./5.;
  char energy[10] = "";
   
  /* other values */
  int i, j, k, l;
  
  int seed = -1; /*negative integer */
  
  double *x, *y, *z, *xall, *Gamma, *H, *J;
  char matfilename[50], numberstr[8];
  
  strcpy(energy, "neutral"); 
  for (i=1; i<argc; ++i) {
    if (!strcmp(argv[i],"-help") || !strcmp(argv[i],"-h")) {
      printf ("usage:  %s [options] \n"
              "where options are:\n"
              "-energy    1 lowest, 4 neutral, 7 highest\n"
              "-M         number of strong vofrtices\n"
              "-N         number of total vofrtices\n"
              ,argv[0]
      );
      exit (1);
    }
    if (!strcmp(argv[i],"-energy")) {
      strcpy(energy,argv[++i]);
    }
    if (!strcmp(argv[i],"-M")) {
      M = atoi(argv[++i]);
      printf("Using %d strong vortices\n", M);
    }
    if (!strcmp(argv[i],"-N")) {
      N = atoi(argv[++i]);
      printf("Using %d total vortices\n", N);
    }
  }
    
  if(!strcmp(energy,"lowest")) {
    ENtarget = -2.;
  }
  else if(!strcmp(energy,"lower")) {
    ENtarget = -1.;
  }
  else if(!strcmp(energy,"low")) {
    ENtarget = -.5;
  }
  else if(!strcmp(energy,"neutral")) {
    ENtarget = 0.;
  }
  else if(!strcmp(energy,"high")) {
    ENtarget = .5;
  }
  else if(!strcmp(energy,"higher")) {
    ENtarget = 1.;
  }
  else if(!strcmp(energy,"highest")) {
    ENtarget = 2.;
  }
  else {
    printf("unknown energy specified, reverting to \"neutral\n");
    strcpy(energy, "neutral");
    ENtarget = 0.;
  }

  /* memory allocation */
  x     = (double *) malloc(N*(sizeof(double)));
  y     = (double *) malloc(N*(sizeof(double)));
  z     = (double *) malloc(N*(sizeof(double)));
  xall  = (double *) malloc(3*N*(sizeof(double)));
  Gamma = (double *) malloc(N*(sizeof(double)));
  H     = (double *) malloc(1*(sizeof(double)));
  J     = (double *) malloc(3*(sizeof(double)));
  
#if GAMMA == 1
  for (i=0; i<M/2; ++i) {
    Gamma[i]     =  GammaA;
    Gamma[M/2+i] = -GammaA;
  }
  for (i=0; i<(N-M)/2; ++i) {
    Gamma[M+i]       =  GammaB;
    Gamma[(N+M)/2+i] = -GammaB;
  }
#elif GAMMA == 2
  skewGamma(Gamma,x,y, lambda, N, &seed);
#elif GAMMA == 3
  for (i=0; i<N; ++i) {
    Gamma[i]     =  GammaA;
  }
#else

#endif
  //readTargets(&ENtarget, &Jtarget, energy);
  printf("Generating initial condition(s) with %d vortices with energy %f and momentum (0,0,%f)\n", N, ENtarget, Jtarget);
  for (l=L0; l<L; ++l) 
  {
    seed -= l;
    
    sprintf(numberstr,"%d",l);
    //sprintf(matfilename, "ic/Lagrange/%s", energy);
    sprintf(matfilename, "../ic/ics_%s_M%dN%d", energy, M, N);
    if(L>1)
      strncat(matfilename, numberstr,strlen(numberstr));
    strncat(matfilename, ".mat",4); 
    fprintf(stderr,"MatLab file will be opened (later) at %s\n", matfilename);  
  
        
    /* initial placement */
    initialPlace(x, y, z, 0, N, Jtarget, Gamma, &seed);
    printf("initial places placed (%f,%f,%f) (%f,%f)\n", x[0], y[0], z[0], Gamma[0], Gamma[N-1]);

    for(i=0;i<2;++i)
      perturbPlace(x, y, z, Gamma, 1.0, H, 0, N, &seed);
    printf("initial places placed and perturbed (%f,%f,%f) (%f,%f)\n", x[0], y[0], z[0], Gamma[0], Gamma[N-1]);
    
#if GAMMA > 1
    while(fabs(J[0])+fabs(J[1])+fabs(J[2]-Jtarget) > sqrt(tol))
    {
      JiterPlace(x, y, z, Gamma, H, J, Jtarget, 0, N, &seed);
    }
    printf("done with J\n");
#endif
    
#if TARGETMN
    conserved(x, y, z, Gamma, H, J, N);
    printf("Conserved quantities (N): H=%f, J=(%f,%f,%f)\n", H[0], J[0], J[1], J[2]);
    i = 0;
    while(fabs(H[0]-ENtarget) > sqrt(tol))
    {
      HiterPlace(x, y, z, Gamma, H, J, ENtarget, 0, N, &seed);
      ++i;
    }
    /*
    HPlace(x, y, z, Gamma, H, J, ENtarget, tol, 0, N, &seed);
    */
#else
#if TARGETM
    conserved(x, y, z, Gamma, H, J, M);
    printf("Conserved quantities (M): H=%f, J=(%f,%f,%f)\\n", H[0], J[0], J[1], J[2]);
    while(fabs(H[0]-EMtarget) > sqrt(tol))
    {
      HiterPlace(x, y, z, Gamma, H, J, EMtarget, 0, M, &seed);
    }
#endif
    
    conserved(x, y, z, Gamma, H, J, N);
    printf("Conserved quantities (N): H=%f, J=(%f,%f,%f)\\n", H[0], J[0], J[1], J[2]);
#if TARGETN
    while(fabs(H[0]-ENtarget) > tol)
    {
      HiterPlace(x, y, z, Gamma, H, J, ENtarget, M, N, &seed);
    }
#else
    for(i=0; i<10; ++i)
    {
      perturbPlace(x, y, z, Gamma, 1.0, H, 0, N, &seed);
    }
#endif
#endif

    /* storing etc. */
    for(i=0;i<N;++i) {
      xall[3*i]   = x[i];
      xall[3*i+1] = y[i];
      xall[3*i+2] = z[i];
    }
    
    MATFile *pmatfile = matOpen(matfilename,"w"); 
    writeArrayToMatlab(pmatfile, xall,  "X0",    3, N);
    writeArrayToMatlab(pmatfile, Gamma, "Gamma", N, 1);
    conserved(x, y, z, Gamma, H, J, M);
    printf("Conserved quantities (M): H=%e, J=(%e,%e,%e)\n", H[0], J[0], J[1], J[2]);
    writeArrayToMatlab(pmatfile, H,  "HM0",      1, 1);
    writeArrayToMatlab(pmatfile, J,  "JM0",      3, 1);
    conserved(x, y, z, Gamma, H, J, N);
    printf("Conserved quantities (N): H=%e, J=(%e,%e,%e)\n", H[0], J[0], J[1], J[2]);
    writeArrayToMatlab(pmatfile, H,  "H0",       1, 1);
    writeArrayToMatlab(pmatfile, J,  "J0",       3, 1);    
    matClose(pmatfile); /* closes the .mat file */
    fprintf(stderr,"MatLab file stored at %s\n\n", matfilename);
  }
}

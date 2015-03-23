/* computes the vortices on a sphere using a barrier-parallelism and higher-order composite methods */
// #define _XOPEN_SOURCE 600

#define _REENTRANT
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mat.h>
#include "normal.h"
#include <dispatch/dispatch.h>
#include <sys/time.h>

//#define THREADS 3
#define HAMILTONIAN 1
#define OMEGA 1
#define AMATRIX 0
#define DEBUG 0
#define MATLABICS 0
#define TARGETICS 0
#define STOREALLX 0
#define STOREALLH 0
#define HFAC 1

void conserved(double *x, double *y, double *z, double *Gamma, double *H, double *J, int thread);

int Ngroup;
int N, M=4, Nstore, Nsteps, Nstages, order = 2; /* 7 is a large 6 */
int interval = 100;
double Tend = 1, dt0 = 0.0001, Etarget = 10.;
char *matfilename = "time.mat";

dispatch_semaphore_t  sema[  THREADS];
dispatch_semaphore_t gsema[2*THREADS*(4*THREADS-2)];

#if AMATRIX
int *A_g_test;
#endif

typedef struct {
  char name[100];
  int size;
  double lower;
  double upper;
  double invdelta;
  double *histogram;
} hist;

int weights(double *w, int order);

void writeArrayToMatlab(MATFile *pmatfile, double *store, char *name, int m, int n)
{
  mxArray *matlabArray;
  
  matlabArray = (mxArray *) mxCreateDoubleMatrix(m,n,mxREAL); /* creates a mxn MATLAB array of reals */
  memcpy((void *)(mxGetPr(matlabArray)),   (void *)store, m*n*sizeof(double)); /* copies the memory (m*n) at Q to tseries */
  matPutVariable(pmatfile,name,matlabArray); /* puts the array matlabx into the .mat file as "dt" */
  mxDestroyArray(matlabArray); /* clears the array (for memory I presume) */  
  printf("Written %s to Matlabfile, it was %d by %d\n", name, m, n);
}

void randPlace(double *x, double *y, double *z, int M, int N, int *pidum)
{
  int i, j=0;
  double *theta, pi=acos(-1);
  theta = r8vec_uniform_01_new(2*(N-M), pidum);
  
  for (i=M; i<N; ++i) {
    theta[j]   = 2*pi*theta[j];    
    theta[j+1] = acos(1-2*theta[j+1]);
    x[i] = sin(theta[j])*sin(theta[j+1]);
    y[i] = cos(theta[j])*sin(theta[j+1]);
    z[i] =               cos(theta[j+1]);
    j   += 2;
  }
}

void targetPlace(double *x, double *y, double *z, double *Gamma, int M, int N, int *pidum)
{
  int i, j, l;
  double E, Eerr = 1.e10, Etemp;
  double *u, *v, *w, pi = acos(-1);
  
  u = (double *) malloc(N*(sizeof(double)));
  v = (double *) malloc(N*(sizeof(double)));
  w = (double *) malloc(N*(sizeof(double)));
  
  for(i=0; i<M; ++i) {
    u[i] = x[i];
    v[i] = y[i];
    w[i] = z[i];
  }
    
  randPlace(x, y, z, M, N, pidum);
  
  for(l=0; l<N; ++l) {
    E = 0;
    randPlace(u, v, w, M, N, pidum);
    for(i=1; i<N; ++i) {
      for(j=0; j<i; ++j) {
        E -= Gamma[i]*Gamma[j]*log( 2 * ( 1. - u[i]*u[j] - v[i]*v[j] - w[i]*w[j] ) );
      }
    }
    E /= 4*pi;
    Etemp = (E-Etarget > 0) ? E-Etarget : Etarget-E;
    if(Etemp < Eerr) {
      Eerr = Etemp;
      for(i=0; i<N; ++i) {
        x[i] = u[i];
        y[i] = v[i];
        z[i] = w[i];
      }
    }
  }
  
  E = 0;
  for(i=1; i<N; ++i) {
    for(j=0; j<i; ++j) {
      E -= Gamma[i]*Gamma[j]*log( 2 * ( 1. - x[i]*x[j] - y[i]*y[j] - z[i]*z[j] ) );
    }
  }
  printf("E: %f\n", E);
  
}



void barrier_wait_all(int thread)
{
  int i;
  
#if DEBUG > 1
  printf("Thread %d waiting!\n", thread);
#endif
  for(i=0;i<THREADS;++i) {
    if(thread != i) {
      /* signal all other threads "I'm ready"*/
      dispatch_semaphore_signal(sema[i]);
    }
  }
  for(i=0;i<THREADS;++i) {
    if(thread != i) {
      /* wait for all other threads to be ready*/
      dispatch_semaphore_wait(sema[thread], DISPATCH_TIME_FOREVER );
    }
  }
}

void barrier_signal(int i, int j, int k, int thread)
{
#if DEBUG > 1
  printf("Thread %d signaling %d and %d are available for round %d. (%d, %d)\n", thread, i, j, k, 2*THREADS*k + i, 2*THREADS*k + j);
#endif
  /* signal all other threads "I'm ready with i and j"*/
  dispatch_semaphore_signal(gsema[2*THREADS*k + i]);
  dispatch_semaphore_signal(gsema[2*THREADS*k + j]);
#if DEBUG > 10
  printf("Thread %d signaled %d and %d are available. (%d, %d)\n", thread, i, j, 2*THREADS*k + i, 2*THREADS*k + j);
#endif
}

void barrier_wait(int i, int j, int k, int thread)
{
#if DEBUG > 1
  printf("Thread %d waiting for %d and %d at level %d. (%d, %d)\n", thread, i, j, k, 2*THREADS*k + i, 2*THREADS*k + j);
#endif
  /* wait for i and j to become available */
  dispatch_semaphore_wait(gsema[2*THREADS*k + i], DISPATCH_TIME_FOREVER );
  dispatch_semaphore_wait(gsema[2*THREADS*k + j], DISPATCH_TIME_FOREVER );
#if DEBUG > 10
  printf("Thread %d waited for %d and %d! (%d, %d)\n", thread, i, j, 2*THREADS*k + i, 2*THREADS*k + j);
#endif
}

void barrier_signal_single(int i, int k, int thread)
{
#if DEBUG > 1
  printf("Thread %d signaling %d is available for round %d. (%d)\n", thread, i, k, 2*THREADS*k + i);
#endif
  /* signal all other threads "I'm ready with i"*/
  dispatch_semaphore_signal(gsema[2*THREADS*k + i]);
#if DEBUG > 10
  printf("Thread %d signaled %d is available. (%d)\n", thread, i, 2*THREADS*k + i);
#endif
}

void barrier_wait_single(int i, int k, int thread)
{
#if DEBUG > 1
  printf("Thread %d waiting for %d at level %d. (%d)\n", thread, i, k, 2*THREADS*k + i);
#endif
  /* wait for i to become available */
  dispatch_semaphore_wait(gsema[2*THREADS*k + i], DISPATCH_TIME_FOREVER );
#if DEBUG > 10
  printf("Thread %d waited for %d! (%d)\n", thread, i, 2*THREADS*k + i);
#endif
}



void lonestep(double *x, double *y, double *z, double *Gamma, double dt, int gi)
{
  /* Rotation of the sphere */
  double s, c, xtemp;
  int i;
  
#if OMEGA
  for(i=gi*Ngroup;i<(gi+1)*Ngroup;++i) {
    s = sin(dt*OMEGA);
    c = cos(dt*OMEGA);
    xtemp = x[i]*c - y[i]*s;
    y[i]  = x[i]*s + y[i]*c;
    x[i]  = xtemp;
  }
#endif
}

void pairstep(double *x, double *y, double *z, double *Gamma, double dt, int i, int j)
{
  double cx, cy, cz, cnorm, s, c, cc1, cc2, xtemp, ytemp;
  
#if HFAC
  double pi= acos(-1);
  double Hloc = 4*pi* (1.0 - x[i]*x[j] - y[i]*y[j] - z[i]*z[j]);
#else
  double Hloc = 1.0 - x[i]*x[j] - y[i]*y[j] - z[i]*z[j];
#endif
  
#if DEBUG > 2
  printf("%3d %3d (%f,%f,%f) (%f,%f,%f)\n", i, j, x[i], y[i], z[i], x[j], y[j], z[j]);
#endif
#if AMATRIX
  A_g_test[N*i+j] += 1;
  A_g_test[N*j+i] += 1;
#endif
  
  cx = 0.5*dt* (Gamma[i]*x[i] + Gamma[j]*x[j])/Hloc;
  cy = 0.5*dt* (Gamma[i]*y[i] + Gamma[j]*y[j])/Hloc;
  cz = 0.5*dt* (Gamma[i]*z[i] + Gamma[j]*z[j])/Hloc;
  
  cnorm = sqrt(cx*cx + cy*cy + cz*cz);
  
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
  
  /*Hloc -= 1.0 - x[i]*x[j] - y[i]*y[j] - z[i]*z[j];
  printf("%e\n", Hloc);*/  
#if DEBUG > 2
  printf("%3d %3d (%f,%f,%f) (%f,%f,%f) cnorm: %f (%f,%f,%f) (%f,%f)\n", i, j, x[i], y[i], z[i], x[j], y[j], z[j], cnorm, cx, cy, cz, Gamma[i], Gamma[j]);
#endif
}

void groupstep(double *x, double *y, double *z, int gi, int gj, double *Gamma, double dt, int thread)
{
  int i, j, k, l;

#if DEBUG
  printf("Groupstep for groups %3d and %3d on thread %d\n", gi, gj, thread);
#endif
  
  /* Round robin WITHIN two separate! groups way up */
  for (k=0;k<Ngroup-1;++k)
  {
    i = Ngroup-1;
    j =        k;
    pairstep(x, y, z, Gamma, dt, i+Ngroup*gi, j+Ngroup*gi);
    pairstep(x, y, z, Gamma, dt, i+Ngroup*gj, j+Ngroup*gj);
    for (l=1;l<Ngroup/2;++l) {
      i = (Ngroup-1-l+k)%(Ngroup-1);
      j =          (l+k)%(Ngroup-1);
      pairstep(x, y, z, Gamma, dt, i+Ngroup*gi, j+Ngroup*gi);
      pairstep(x, y, z, Gamma, dt, i+Ngroup*gj, j+Ngroup*gj);
    }
  }
  
  /* stuff on a single vortex */
  lonestep(x, y, z, Gamma, dt, gi);
  lonestep(x, y, z, Gamma, dt, gj);
    
  
  /* Round robin WITHIN two separate! groups way down */
  for (k=Ngroup-2;k>=0;--k)
  {
    i = Ngroup-1;
    j =        k;
    pairstep(x, y, z, Gamma, dt, i+Ngroup*gi, j+Ngroup*gi);
    pairstep(x, y, z, Gamma, dt, i+Ngroup*gj, j+Ngroup*gj);
    for (l=1;l<Ngroup/2;++l) {
      i = (Ngroup-1-l+k)%(Ngroup-1);
      j =          (l+k)%(Ngroup-1);
      pairstep(x, y, z, Gamma, dt, i+Ngroup*gi, j+Ngroup*gi);
      pairstep(x, y, z, Gamma, dt, i+Ngroup*gj, j+Ngroup*gj);
    }
  }
}

void battlestep(double *x, double *y, double *z, double *Gamma, double dt, int gi, int gj, int dir, int thread)
{
  int l = thread;
  int i, j, m, n;
#if DEBUG
  printf("Battling groups %3d and %3d on thread %d\n", gi, gj, thread);
#endif

  if(dir > 0)
  {
    for(m=0;m<Ngroup;m++)
    {
      for(n=0;n<Ngroup;n++)
      {
        i =    n        +Ngroup*gi;
        j = (m+n)%Ngroup+Ngroup*gj;
        pairstep(x, y, z, Gamma, dt, i, j);
      }
    }
  } else {
    for(m=Ngroup-1;m>=0;m--)
    {
      for(n=0;n<Ngroup;n++)
      {
        i =    n        +Ngroup*gi;
        j = (m+n)%Ngroup+Ngroup*gj;
        pairstep(x, y, z, Gamma, dt, i, j);
      }
    }
  }
}

#if 0
void roundstep(double *x, double *y, double *z, double *Gamma, double dt, int k, int dir, int thread)
{
  int l = thread, gi, gj;
  
  if(l == 0) {
    gi = 2*THREADS-1;
    gj = k;
  } else {
    gi = (2*THREADS-1-l+k)%(2*THREADS-1);
    gj =             (l+k)%(2*THREADS-1);
  }
  
  barrier_wait(gi, gj, (k*dir + 4*THREADS-2 + (dir-1)/2 )%(4*THREADS-2), thread); /* wait for prev step */
  battlestep(x, y, z, Gamma, dt, gi, gj, dir, thread);
  barrier_signal(gi, gj, (k*dir + 4*THREADS-1 + (dir-1)/2 )%(4*THREADS-2), thread); /* signal for next step */
}
#endif

void symmetricstep(double *x, double *y, double *z, double *Gamma, double dt, int thread)
{
  int gi, gj, k, l = thread, dir;
  
  if(thread == 0) {
    /* on first thread!!! */
    gi = 2*THREADS-1;
   
    /* Round robin OF/BETWEEN groups (way up) */
    dir = 1;
    barrier_wait_single(gi, 0, thread);
    for (gj=0; gj<2*THREADS-2; ++gj) {
      barrier_wait_single(gj, gj, thread);
      battlestep(x, y, z, Gamma, dt, gi, gj, dir, thread);
      barrier_signal_single(gj, gj+1, thread);
    }
    gj = 2*THREADS-2;
    barrier_wait_single(gj, gj, thread);
    battlestep(x, y, z, Gamma, dt, gi, gj, dir, thread);

    /* computing the interaction of vortices from a single group */
    groupstep(x, y, z, gi, gj, Gamma, dt, thread);
    
    /* Round robin OF/BETWEEN groups (way down) */
    dir = -1;
    battlestep(x, y, z, Gamma, dt, gi, gj, dir, thread);
    for (gj=2*THREADS-3; gj>=0; --gj) {
      barrier_signal_single(gj+1, 4*THREADS-4-gj, thread);
      barrier_wait_single(gj, 4*THREADS-4-gj, thread);
      battlestep(x, y, z, Gamma, dt, gi, gj, dir, thread);
    }
    barrier_signal(gi, 0, 0, thread);
    
  } else {
    /* on other threads */
    gi = l;
    gj = 2*THREADS-1-l;
   
    /* Round robin OF/BETWEEN groups (way up) */
    dir = 1;
    barrier_wait(gi, gj, 0, thread);
    for (k=0; k<2*THREADS-2; ++k) {
      battlestep(x, y, z, Gamma, dt, gi, gj, dir, thread);
      if((k+1-l)%(THREADS-1) == 0) {
        barrier_signal(gi, gj, k+1, thread);
        gi = (gi+THREADS-1)%(2*THREADS-1);
        gj = (gi+1)%(2*THREADS-1);
        barrier_wait(gi, gj, k+1, thread);
      } else {
        barrier_signal_single(gj, k+1, thread);
        gj = (gj+2)%(2*THREADS-1);
        barrier_wait_single(gj, k+1, thread);
      }
    }
    battlestep(x, y, z, Gamma, dt, gi, gj, dir, thread);
    
    /* computing the interaction of vortices from a single group */
    groupstep(x, y, z, gi, gj, Gamma, dt, thread);
    
    /* Round robin OF/BETWEEN groups (way down) */
    dir = -1;
    battlestep(x, y, z, Gamma, dt, gi, gj, dir, thread);
    for (k=2*THREADS-3; k>=0; --k) {
      if((k+1-l)%(THREADS-1) == 0) {
        barrier_signal(gi, gj, 4*THREADS-4-k, thread);
        gi = (gi+THREADS)%(2*THREADS-1);
        gj = (gi+2*THREADS-3)%(2*THREADS-1); /* actually -2, but to avoid negative ints... */
        barrier_wait(gi, gj, 4*THREADS-4-k, thread);
      } else {
        barrier_signal_single(gj, 4*THREADS-4-k, thread);
        gj = (gj+2*THREADS-3)%(2*THREADS-1); /* actually -2, but to avoid negative ints... */
        barrier_wait_single(gj, 4*THREADS-4-k, thread);
      }
      battlestep(x, y, z, Gamma, dt, gi, gj, dir, thread);
    }
    barrier_signal(gi, gj, 0, thread);
  }
}

void spherestep(double *x, double *y, double *z, double *Gamma, double *w, int Nstages, int thread)
{
  int i;
  double dt;
  
  for(i=0; i<Nstages; ++i) {
    dt = dt0 * w[i];
    symmetricstep(x, y, z, Gamma, dt, thread);
  }
}

void conserved(double *x, double *y, double *z, double *Gamma, double *H, double *J, int thread)
{
#if DEBUG
  printf("Computing the conserved quantities\n");
#endif
  int l = thread;
  int i, j, k;
  
  H[l]   = 0.;
  J[3*l]   = 0.;
  J[3*l+1] = 0.;
  J[3*l+2] = 0.;
  
  for(i=l*Ngroup; i<(l+1)*Ngroup; ++i)
  {
    /* do one "i" from the beginning, a small one */
#if 0
    H[l+1] += Gamma[i] * OMEGA * z[i];
#endif
    for(j=0; j<i; ++j) {
      H[l]   -= Gamma[i]*Gamma[j]*log( 2 * ( 1. - x[i]*x[j] - y[i]*y[j] - z[i]*z[j] ) );
    }
    J[3*l]   += Gamma[i]*x[i];
    J[3*l+1] += Gamma[i]*y[i];
    J[3*l+2] += Gamma[i]*z[i];
    
    /* do one "i" from the end, a large one */
    k = N - i - 1;
#if 0
    H[l+1] += Gamma[k] * OMEGA * z[k];
#endif
    for(j=0; j<k; ++j) {
      H[l] -= Gamma[k]*Gamma[j]*log( 2 * ( 1. - x[k]*x[j] - y[k]*y[j] - z[k]*z[j] ) );
    }
    J[3*l]   += Gamma[k]*x[k];
    J[3*l+1] += Gamma[k]*y[k];
    J[3*l+2] += Gamma[k]*z[k];
  }
}

void energy_str(double *x, double *y, double *z, double *Gamma, double *HM, double *HI, double *JM)
{
  double pi = acos(-1);
  int i, j, k;
  
  for(i=0; i<M; ++i)
  {
    JM[0] += Gamma[i]*x[i];
    JM[1] += Gamma[i]*y[i];
    JM[2] += Gamma[i]*z[i];
    for(j=0; j<i; ++j) {
      HM[0] -= Gamma[i]*Gamma[j]*log( 2 * ( 1. - x[i]*x[j] - y[i]*y[j] - z[i]*z[j] ) );
    }
    for(j=M; j<N; ++j) {
      HI[0] -= Gamma[i]*Gamma[j]*log( 2 * ( 1. - x[i]*x[j] - y[i]*y[j] - z[i]*z[j] ) );
    }
  }
#if HFAC
  HM[0] /= 4*pi;
  HI[0] /= 4*pi;
#endif
}

#if STOREALLH
void output(double *x, double *y, double *z, double *Gamma, double *H, double *J, double *xstore, double *Hstore, double *Hpstore, double *Jstore, int thread, int stepStore)
#else
void output(double *x, double *y, double *z, double *H, double *J, double *xstore, double *Hstore, double *Jstore, int thread, int stepStore)
#endif
{
  double pi = acos(-1);
  int l,m, k = thread, Np = N*N;
#if STOREALLX > 0
  for (l=2*k*Ngroup; l<2*(k+1)*Ngroup; ++l) {
    xstore[3*N*stepStore+l]     = x[l];
    xstore[3*N*stepStore+l+N]   = y[l];
    xstore[3*N*stepStore+l+2*N] = z[l];
  }
#endif
#if STOREALLH
  for (l=k*Ngroup; l<(k+1)*Ngroup; ++l) {
    for (m=0; m<l; ++m) {
      Hpstore[Np*stepStore+l*N+m] = -Gamma[l]*Gamma[m]*log( 2 * ( 1. - x[l]*x[m] - y[l]*y[m] - z[l]*z[m] ) )/(4*pi);
    }
  }
  for (l=(2*THREADS-k-1)*Ngroup; l<(2*THREADS-k)*Ngroup; ++l) {
    for (m=0; m<l; ++m) {
      Hpstore[Np*stepStore+l*N+m] = -Gamma[l]*Gamma[m]*log( 2 * ( 1. - x[l]*x[m] - y[l]*y[m] - z[l]*z[m] ) )/(4*pi);
    }
  }
#endif
  
  if(thread == 0) {
#if STOREALLX < 1
    for (l=0; l<M; ++l) {
      xstore[3*M*stepStore+l]     = x[l];
      xstore[3*M*stepStore+l+M]   = y[l];
      xstore[3*M*stepStore+l+2*M] = z[l];
    }
#endif
    for(l=0;l<THREADS;++l) {
      Jstore[3*stepStore]   += J[3*l]; 
      Jstore[3*stepStore+1] += J[3*l+1]; 
      Jstore[3*stepStore+2] += J[3*l+2]; 
      Hstore[stepStore]     += H[l];
    }
#if HFAC
  Hstore[stepStore] /= 4*pi;
#endif
  }
}

static void parallel (void *workspace, size_t thread) {
  double *x       = workspace;
  double *y       = x + N;
  double *z       = y + N;
  double *Gamma   = z + N;
  double *H       = Gamma + N;
  double *J       = H + THREADS;
  double *xstore  = J + 3*THREADS;
  double *Hstore  = xstore + 3*N*Nstore;
  double *HMstore = Hstore + Nstore;
  double *HIstore = HMstore + Nstore;
  double *Jstore  = HIstore + Nstore;
  double *JMstore = Jstore  + 3*Nstore;
  double *w       = JMstore + 3*Nstore;
#if STOREALLH
  double *Hpstore = w + 100;
#endif
  
  double simutime = 0.;
  
  int i, step, stepStore = 0, Nperthread = Ngroup;
  
  /* signal initially that all threads are ready */
  barrier_signal(2*thread, 2*thread+1, 0, thread);
  
  /* main code running in parallel! */
  for(step=0; step<Nsteps; ++step)
  {
    spherestep(x, y, z, Gamma, w, Nstages, thread);
    simutime += dt0;
    if( step%interval == interval-1)
    {
      barrier_wait_all(thread);
      conserved(x, y, z, Gamma, H, J, thread);
      barrier_wait_all(thread);
#if STOREALLH
      output(x, y, z, Gamma, H, J, xstore, Hstore, Hpstore, Jstore, thread, stepStore);
#else
      output(x, y, z, H, J, xstore, Hstore, Jstore, thread, stepStore);
#endif
      
      if(thread == 0) {
        /*printf("Step %d of %d. Simutime: %f, H: %f\n", step, Nsteps, simutime, Hstore[stepStore]);*/
        printf("Step %d of %d. Simutime: %3f, H: %5f, J: (%5f,%5f,%5f) Jxy: %f\n", step, Nsteps, simutime, Hstore[stepStore], Jstore[3*stepStore], Jstore[3*stepStore+1], Jstore[3*stepStore+2], Jstore[3*stepStore]*Jstore[3*stepStore] + Jstore[3*stepStore+1]*Jstore[3*stepStore+1] );
        energy_str(x, y, z, Gamma, HMstore+stepStore, HIstore+stepStore, JMstore+3*stepStore);
      }
      ++stepStore;
      Hstore[stepStore] = 0.;
    }
    barrier_wait_all(thread);
  }
  
}

int main (int argc, char *argv[])
{  
  
  int i, j;
  int idum = -1; /* negative integer */
  
  for(i=0;i<THREADS;++i) {
    sema[i] = dispatch_semaphore_create(0);
  }
  printf("Nsteps: %d\n", Nsteps);
  for(i=0;i<2*THREADS*(4*THREADS-1);++i) {
    gsema[i] = dispatch_semaphore_create(0);
  }
  printf("Nsteps: %d\n", Nsteps);
  
  for (i=1; i<argc; ++i) {
    if (!strcmp(argv[i],"-help") || !strcmp(argv[i],"-h")) {
      printf ("usage:  %s [options] \n"
              "where options are:\n"
              "-g #     number of vortices in half a group\n"
              ,argv[0]
    );
      exit (1);
    }
    if (!strcmp(argv[i],"-g")) Ngroup = 2 * atoi(argv[++i]);
  }
  
  dispatch_queue_t queue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0);

  struct timeval  tv1, tv2;
  double time;
  
  N = Ngroup * 2 * THREADS;
  Nsteps = ceil(Tend/dt0);
  Nstore = ceil(Tend/(interval*dt0));
           
  int Nwork = 4*N + 4*THREADS + (3 + 3*N + 2*3) * Nstore + 100; /* pos,G + H,J + (Hs,xs,Js) + weights */
#if STOREALLH
  Nwork += N*N * Nstore;
#endif
  
  double *workspace;
  workspace = (double *) malloc(Nwork*(sizeof(double)));
  double *x       = workspace;
  double *y       = x + N;
  double *z       = y + N;
  double *Gamma   = z + N;
  double *H       = Gamma + N;
  double *J       = H + THREADS;
  double *xstore  = J + 3*THREADS;
  double *Hstore  = xstore + 3*N*Nstore;
  double *HMstore = Hstore + Nstore;
  double *HIstore = HMstore + Nstore;
  double *Jstore  = HIstore + Nstore;
  double *JMstore = Jstore  + 3*Nstore;
  double *w       = JMstore + 3*Nstore;
#if STOREALLH
  double *Hpstore = w + 100;
#endif
  
  Nstages = weights(w, order);
  
#if MATLABICS
  readInitCondFromMatLab(x, y, z, Gamma, N, "../ics_highest.mat");
#else
  for (i=0; i<M/2; ++i) {
    Gamma[i]     =  5.;
    Gamma[M/2+i] = -5.;
  }
  for (i=0; i<(N-M)/2; ++i) {
    Gamma[M+i]       =  1.;
    Gamma[(N+M)/2+i] = -1.;
  }
  
  x[0] =-1.0;
  x[1] = 1.0;
  y[2] = 1.0; 
  y[3] =-1.0;
  y[0] = y[1] = x[2] = x[3] = 0.0;
  z[0] = z[1] = z[2] = z[3] = 0.0;
#if TARGETICS
  targetPlace(x, y, z, Gamma, M, N, &idum);
#else
  randPlace(x, y, z, M, N, &idum);
#endif
#endif
  
#if AMATRIX
  A_g_test = (int *) malloc(N*N*(sizeof(int)));
  for(i=0; i<N; ++i)
    for(j=0; j<N; ++j)
      A_g_test[i*N+j] = 0;
#endif
  printf("Started for N = %d, Nsteps = %d.\n", N, Nsteps);
  fprintf(stderr,"MatLab file will be opened (later) at %s\n", matfilename);  
  
  gettimeofday(&tv1, NULL);    
  dispatch_apply_f(THREADS, queue, workspace, parallel);
  gettimeofday(&tv2, NULL);
  
  time = (double) (tv2.tv_usec - tv1.tv_usec)/1000000 + (double) (tv2.tv_sec - tv1.tv_sec);
  printf ("\nTotal time = %f seconds\n\n", time);
  
  FILE *file = fopen("onefile.txt", "a");
  fprintf(file, "%d, %d, %f, %d, %d, %f\n", THREADS, Ngroup, time, Nsteps, interval, dt0);
  fclose(file);
  
#if AMATRIX
  for(i=0; i<N; ++i)
  {
    printf("[");
    for(j=0; j<N; ++j)
    {
      printf("%5d", A_g_test[i*N+j]);
    }
    printf("]\n");
  }
#endif

  MATFile *pmatfile = matOpen(matfilename,"w"); 
      
  writeArrayToMatlab(pmatfile, Gamma,  "G", N,   1);
#if STOREALLX
  writeArrayToMatlab(pmatfile, xstore, "x", N*3, Nstore);
#else
  writeArrayToMatlab(pmatfile, xstore, "x", M*3, Nstore);
#endif
  writeArrayToMatlab(pmatfile, Hstore, "H", 1,   Nstore);
#if STOREALLH
  writeArrayToMatlab(pmatfile, Hpstore,"Hp",N*N, Nstore);
#endif
  writeArrayToMatlab(pmatfile, HMstore,"HM",1,   Nstore);
  writeArrayToMatlab(pmatfile, HIstore,"HI",1,   Nstore);
  writeArrayToMatlab(pmatfile, Jstore, "J", 3,   Nstore);
  writeArrayToMatlab(pmatfile, JMstore,"JM",3,   Nstore);
  writeArrayToMatlab(pmatfile, &dt0,   "dt",1,   1);
  fprintf(stderr,"MatLab file stored at %s\n", matfilename);
  
  matClose(pmatfile); /* closes the .mat file */
  
  printf("Finished for N = %d, Nsteps = %d.\n", N, Nsteps);
}

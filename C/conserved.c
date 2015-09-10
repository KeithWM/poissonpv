#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mat.h>
#include "normal.h"
//#include <dispatch/dispatch.h>
#include <sys/time.h>
#include "parallel.h"

double Hpair(int i, int j, struct Solution *psol, struct Parameters parc, struct Dynamics dync) {
  double pi = acos(-1);

  if(i < parc.N && j < parc.N) {
//    printf("[%d,%d]\n", i,j);
    return -dync.G[i]*dync.G[j]*log( 2 * ( 1. - psol->x[i]*psol->x[j] - psol->y[i]*psol->y[j] - psol->z[i]*psol->z[j] ) )/(4*pi);
  } else {
    return 0.;
  }
}

void incHbetween(double *Hct, unsigned char k, unsigned char l, struct Solution *psol, struct Parameters parc, struct Dynamics dync) {
  int i,j,c = 0;
  
//  printf("%d: [%d,%d]\n",c, k, l);
  while((c = pairbetween(&i, &j, c, parc.L, 1)) <= parc.L*parc.L) {
    *Hct+= Hpair(k*parc.L+i, l*parc.L+j, psol, parc, dync);
  }
//   c-= 2;
//   printf("c: %d\n", c);
//   while((c = pairbetween(&i, &j, c, parc.L, -1)) >= -1) {
//     printf("%d: [%d,%d]\n",c, k*parc.L+i, l*parc.L+j);
//   }
//   printf("c: %d\n", c);
}

void incHwithin(double *Hct, unsigned char k, unsigned char l, struct Solution *psol, struct Parameters parc, struct Dynamics dync) {
  int i,j,c = 0;
  
//   printf("%d: [%d,%d]\n",c, k, l);
  while((c = pairwithin(&i, &j, c, parc.L, 1)) <= parc.L*(2*parc.L-1)) {
    i = (i < parc.L) ? k*parc.L+i : (l-1)*parc.L+i;
    j = (j < parc.L) ? k*parc.L+j : (l-1)*parc.L+j;
    *Hct+= Hpair(i, j, psol, parc, dync);
  }
}


void incJ(double *Jct, int k, struct Solution *psol, struct Parameters parc, struct Dynamics dync) {
  int i,j;
  for(i=0;i<2*parc.L;++i) {
    j = i+2*k*parc.L;
    Jct[0]+= dync.G[j]*psol->x[j];
    Jct[1]+= dync.G[j]*psol->y[j];
    Jct[2]+= dync.G[j]*psol->z[j];
  }
}

void conserved(struct Solution *psol, struct Parameters parc, struct Dynamics dync, size_t thread) {
  int i,j,k, r,s;
  
  /* compute Hamiltonian */
  psol->Hc[thread] = 0.;
  for(r=0; r<2*THREADS-2; ++r) {
    incHbetween(psol->Hc+thread, dync.A[2*dync.R*thread + r], dync.B[2*dync.R*thread + r], psol, parc, dync);
  }
  incHwithin(psol->Hc+thread, dync.A[2*dync.R*thread + r], dync.B[2*dync.R*thread + r], psol, parc, dync);
//  printf("Hc[%zu]: %.2e\n", thread, psol->Hc[thread]);
  
  /* compute Noether momenta */
  k = thread;
  for(j=0; j<4; ++j) {
    psol->Jc[4*thread+j] = 0.;
    incJ(psol->Jc+4*thread, k, psol, parc, dync);
  }
  
  /* combine  results */
/* DO NOT DO THIS HERE, AS IT REQUIRES WAITING
  if(thread == THREADS-1) {
    psol->H[0] = 0.;
    for(i=0; i<THREADS; ++i) {
      psol->H[0]+= psol->Hc[i];
    }
  }
  if(thread == THREADS-1) {
    psol->J[3] = 0;
    for(j=0; j<3; ++j) {
      psol->J[j] = 0.;
      for(i=0; i<THREADS; ++i) {
        psol->J[j]+= psol->Jc[4*i+j];
      }
      psol->J[3]+= psol->J[j]*psol->J[j];
    }
  }
*/
}

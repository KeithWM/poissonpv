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

void pairstep(int i, int j, struct Solution *psol, struct Parameters parc, struct Dynamics dync, double dt)
{
  double cx, cy, cz, cnorm, s, c, cc1, cc2, xtemp, ytemp;
  double pi= acos(-1);
//  printf("(%d,%d)\n", i,j);
  if(i < parc.N && j < parc.N) {
    double Hloc = 4*pi* (1.0 - psol->x[i]*psol->x[j] - psol->y[i]*psol->y[j] - psol->z[i]*psol->z[j]);
  
    cx = dt* (dync.G[i]*psol->x[i] + dync.G[j]*psol->x[j])/Hloc;
    cy = dt* (dync.G[i]*psol->y[i] + dync.G[j]*psol->y[j])/Hloc;
    cz = dt* (dync.G[i]*psol->z[i] + dync.G[j]*psol->z[j])/Hloc;
    
    cnorm = sqrt(cx*cx + cy*cy + cz*cz);
    
    cx = cx/cnorm;
    cy = cy/cnorm;
    cz = cz/cnorm;
    
    s = sin(cnorm);
    c = 1.-cos(cnorm);
    
    cc1 = cx*psol->x[i] + cy*psol->y[i] + cz*psol->z[i];
    cc2 = cx*psol->x[j] + cy*psol->y[j] + cz*psol->z[j];
    
    xtemp      = psol->x[i] + s * ( cy*psol->z[i]-cz*psol->y[i] ) - c * ( psol->x[i]-cc1*cx );
    ytemp      = psol->y[i] + s * ( cz*psol->x[i]-cx*psol->z[i] ) - c * ( psol->y[i]-cc1*cy );
    psol->z[i] = psol->z[i] + s * ( cx*psol->y[i]-cy*psol->x[i] ) - c * ( psol->z[i]-cc1*cz );
    psol->x[i] = xtemp;
    psol->y[i] = ytemp;
    
    xtemp      = psol->x[j] + s * ( cy*psol->z[j]-cz*psol->y[j] ) - c * ( psol->x[j]-cc2*cx );
    ytemp      = psol->y[j] + s * ( cz*psol->x[j]-cx*psol->z[j] ) - c * ( psol->y[j]-cc2*cy );
    psol->z[j] = psol->z[j] + s * ( cx*psol->y[j]-cy*psol->x[j] ) - c * ( psol->z[j]-cc2*cz );
    psol->x[j] = xtemp;
    psol->y[j] = ytemp;
    
//    printf("[%d,%d], (%f,%f,%f), (%f,%f,%f)\n", i,j, psol->x[i], psol->y[i], psol->z[i], psol->x[j], psol->y[j], psol->z[j]); 
  }
}

void stepbetween(unsigned char k, unsigned char l, struct Solution *psol, struct Parameters parc, struct Dynamics dync, int dir, int *ps, size_t thread) {
  int i,j,c;
  double dt;
  
  // only use full split time step if Lie-Trotter
  dt = (strcmp(parc.compo, "LT") == 0) ? parc.dt_split : .5*parc.dt_split;
  if(parc.synch)
    waitfor(k,l, thread, *ps);
  if(dir == 1) {
    c = 0;
    while((c = pairbetween(&i, &j, c, parc.L,  1)) <= parc.L*parc.L) {
      pairstep(k*parc.L+i,l*parc.L+j, psol, parc, dync, dt);
    }
  } else if(dir == -1) {
    c = parc.L*parc.L-1;
    while((c = pairbetween(&i, &j, c, parc.L, -1)) >= -1) {
      pairstep(k*parc.L+i,l*parc.L+j, psol, parc, dync, dt);
    }
  }
  ++(*ps);
  
  if(parc.synch)
    donewith(k,l, thread, *ps%dync.S);
}

void stepwithin(unsigned char k, unsigned char l, struct Solution *psol, struct Parameters parc, struct Dynamics dync, int dir, int *ps, size_t thread) {
  int i,j,c;
  double dt;
  
  // only use full split time step if Lie-Trotter
  dt = (strcmp(parc.compo, "LT") == 0) ? parc.dt_split : .5*parc.dt_split;
  if(parc.synch)
    waitfor(k,l, thread, *ps);
//   printf("%d: [%d,%d]\n",*ps, k, l);
  if(dir == 1) {
    c = 0;
    while((c = pairwithin(&i, &j, c, parc.L,  1)) <= parc.L*(2*parc.L-1)) {
      i = (i < parc.L) ? k*parc.L+i : (l-1)*parc.L+i;
      j = (j < parc.L) ? k*parc.L+j : (l-1)*parc.L+j;
      pairstep(i,j, psol, parc, dync, dt);
    }
  } else if(dir == -1) {
    c = parc.L*(2*parc.L-1)-1;
    while((c = pairwithin(&i, &j, c, parc.L, -1)) >= -1) {
      i = (i < parc.L) ? k*parc.L+i : (l-1)*parc.L+i;
      j = (j < parc.L) ? k*parc.L+j : (l-1)*parc.L+j;
      pairstep(i,j, psol, parc, dync, dt);
    }
  }
  ++(*ps);
  
  if(parc.synch)
    donewith(k,l, thread, *ps%dync.S);
}

void LTstep(struct Solution *psol, struct Parameters parc, struct Dynamics dync, size_t thread) { // A WHOLE Lie-Trotter splitting time step
  int i, r, s=0;
  int *ps = &s;

  for(r=0; r<2*THREADS-2; ++r) {
    stepbetween(dync.A[2*dync.R*thread + r], dync.B[2*dync.R*thread + r], psol, parc, dync,  1, ps, thread);
  }
  stepwithin(dync.A[2*dync.R*thread + r], dync.B[2*dync.R*thread + r], psol, parc, dync,  1, ps, thread);
}

void Sstep(struct Solution *psol, struct Parameters parc, struct Dynamics dync, size_t thread) { // A WHOLE Strang splitting time step
  int i, r, s=0;
  int *ps = &s;

  for(r=0; r<2*THREADS-2; ++r) {
    stepbetween(dync.A[2*dync.R*thread + r], dync.B[2*dync.R*thread + r], psol, parc, dync,  1, ps, thread);
  }
  stepwithin(dync.A[2*dync.R*thread + r], dync.B[2*dync.R*thread + r], psol, parc, dync,  1, ps, thread);
  stepwithin(dync.A[2*dync.R*thread + r], dync.B[2*dync.R*thread + r], psol, parc, dync, -1, ps, thread);
  for(r=2*THREADS-3; r >= 0; --r) {
    stepbetween(dync.A[2*dync.R*thread + r], dync.B[2*dync.R*thread + r], psol, parc, dync, -1, ps, thread);
  }
}

void spherestep(struct Solution *psol, struct Parameters parc, struct Dynamics dync, size_t thread) {
  int i;

  if(strcmp(parc.compo, "LT") == 0) { // see Hairer, Lubich and Wanner for more details
    LTstep(psol, parc, dync, thread);
  } else if (strcmp(parc.compo, "S") == 0) {
    Sstep(psol, parc, dync, thread);
  } else { // more involved compositions, loop over splits
    for (i=0; i<parc.Nsplit; ++i) {
      parc.dt_split = parc.dts_split[i]*parc.dt;
      Sstep(psol, parc, dync, thread);
    }
  }
}

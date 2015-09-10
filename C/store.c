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

void store(struct Solution *psol, struct Storage *psto, struct Parameters parc, size_t thread, int nop1) {
  int i,j;
  
  if(thread == 0) {
    psto->t[nop1] = psol->t;
  }
 
  /* add H and J components for storage */ 
  if(thread == 0) {
    psto->H[nop1] = 0.;
    for(i=0; i<THREADS; ++i) {
      psto->H[nop1]+= psol->Hc[i];
    }
  }
  if(thread == THREADS-1) {
    psto->J[nop1*4 + 3] = 0;
    for(j=0; j<3; ++j) {
      psto->J[nop1*4 + j] = 0.;
      for(i=0; i<THREADS; ++i) {
        psto->J[nop1*4 + j]+= psol->Jc[4*i+j];
      }
      psto->J[nop1*4 + 3]+= psto->J[nop1*4 + j]*psto->J[nop1*4 + j];
    }
  }
 
  for(i=0;i<2*parc.L;++i) {
    j = i+2*thread*parc.L;
    psto->x[nop1*parc.N + j] = psol->x[j];
    psto->y[nop1*parc.N + j] = psol->y[j];
    psto->z[nop1*parc.N + j] = psol->z[j];
  }
  
//  printf("stored part %zu at position %d\n", thread, nop1);
}

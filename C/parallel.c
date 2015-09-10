#define _REENTRANT
#ifndef __has_extension
#define __has_extension(x) 0
#endif
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mat.h>
#include "normal.h"
#include <dispatch/dispatch.h>
#include <sys/time.h>
#include "parallel.h"
#define SHOWORDER 0

dispatch_semaphore_t  sema[THREADS*THREADS];
dispatch_semaphore_t gsema[2*THREADS*(4*THREADS-2)];

int pairbetween(int *i, int *j, int counter, int L, int dir) {
  // loops over all possible pairs (i,j) with i and j from distinct groups
  *i = counter%L;
  *j = (*i + (counter-*i)/L)%L;
  
  counter+= dir;
  return counter;
}

int pairwithin(int *i, int *j, int counter, int L, int dir) {
  // loops over all possible pairs (i,j) including both sets (i,j) from distinct groups or from the same group
  int k = counter%L, l = (counter - k)/L;
  *i = (k+l)%(2*L-1);
  *j = (k == 0) ? 2*L-1 : (2*L-1-k+l)%(2*L-1);
  
  counter+= dir;
  return counter;
}

void waitfor(int k, int l, size_t thread, int s)
{
  int i;
  char space[50];

#if SHOWORDER
  strcpy(space, "");
  for(i=0; i<thread; ++i) sprintf(space, "%s\t\t", space);
  printf("%s?%zu,%d:%d,%d?\n", space, thread, s, k,l); 
#endif
  /* wait for k and l to become available at step s*/
  dispatch_semaphore_wait(gsema[2*THREADS*s + k], DISPATCH_TIME_FOREVER );
  dispatch_semaphore_wait(gsema[2*THREADS*s + l], DISPATCH_TIME_FOREVER );
#if SHOWORDER
  printf("%s>%zu,%d:%d,%d<\n", space, thread, s, k,l); 
#endif
}

void donewith(int k, int l, size_t thread, int s)
{
  int i;
  char space[50];

#if SHOWORDER
  strcpy(space, "");
  for(i=0; i<thread; ++i) sprintf(space, "%s\t\t", space);
  printf("%s<%zu,%d:%d,%d>\n", space, thread, s, k,l); 
#endif
  /* signal all other threads "I'm ready with k and l at step s"*/
  dispatch_semaphore_signal(gsema[2*THREADS*s + k]);
  dispatch_semaphore_signal(gsema[2*THREADS*s + l]);
}

void barrier(size_t thread)
{
  int i;

  for(i=0;i<THREADS;++i) {
      /* signal all other threads "I'm ready"*/
      dispatch_semaphore_signal(sema[THREADS*thread + i]);
  }
  for(i=0;i<THREADS;++i) {
      /* wait for all other threads to be ready*/
      dispatch_semaphore_wait(sema[THREADS*i + thread], DISPATCH_TIME_FOREVER );
  }
}

void parallel(void **Pass, size_t thread) {
  int ni, no; // inner and outer loops
  struct Parameters *ppar = Pass[0];
  struct Solution   *psol = Pass[1];
  struct Dynamics   *pdyn = Pass[2];
  struct Storage    *psto = Pass[3];
  struct Parameters  parc = *ppar; // copy of par for ease of notation
  struct Dynamics    dync = *pdyn; // copy of dyn for ease of notation
  

  conserved( psol, parc, dync, thread);
  barrier(thread);
  store(psol, psto, parc, thread, 0);
  barrier(thread); // wait for storage of IC before starting
  if(parc.synch)
    donewith(2*thread, 2*thread+1, thread, 0);

  if(thread == THREADS-1) {
    printf("Conserved quantities computed to be H: %.2e,\tJ:(%.2e,%.2e,%.2e),\t|J|^2: %.2e (at start)\n", psto->H[0], psto->J[0], psto->J[1], psto->J[2], psto->J[3]);
  }
  
  for(no=0; no<parc.Nouter; ++no) {
    for(ni=0; ni<parc.Ninner; ++ni) {
      if(thread == 0)
        psol->t+= parc.dt;
      spherestep(psol, parc, dync, thread);
    }
    barrier(thread); // wait for everything before computing conserved quantities and storing solution
    conserved(psol, parc, dync, thread);
    barrier(thread); // wait for everything before computing conserved quantities and storing solution
    store(psol, psto, parc, thread, no+1);
    barrier(thread); // wait 
    if(thread == THREADS-1) {
      printf("t: %2.2f; Conserved quantities computed to be H: %.2e,\tJ:(%.2e,%.2e,%.2e),\t|J|^2: %.2e\n", psto->t[no+1], psto->H[no+1], psto->J[no*4+4], psto->J[no*4+5], psto->J[no*4+6], psto->J[no*4+7]);
//printf("x[%d]: (%f,%f,%f)\n", 0, psol->x[0], psol->x[0], psol->z[0]);
    }
  }
}

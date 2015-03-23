#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mat.h>
#include "parallel.h"

int readInitCondFromMatLab(struct Solution *psol, struct Dynamics *pdyn, struct Parameters par, char *file)
{  
  MATFile *pmat;
  const char *name;
  int   ndim;
  int   i;
  mxArray *pa;
  const int *dims;
  double *X;
  
  dims = (int *) malloc(2*(sizeof(int)));
  
  printf("Reading file %s...\n\n", file);

  pmat = matOpen(file, "r");
  if (pmat == NULL) {
    printf("Error reopening file %s\n", file);
    return(1);
  }

  /* Get headers of all variables */
//   printf("\nExamining the header for each variable:\n");
  psol->t = 0;
  while( (pa = matGetNextVariable(pmat, &name)) != NULL )
  {
    ndim = mxGetNumberOfDimensions(pa);
    dims = mxGetDimensions(pa);
    if( strcmp(name, "X0") == 0 ) {
      X = (double *) malloc(dims[0]*dims[1]*(sizeof(double)));
      memcpy((void *)X, (void *)(mxGetPr(pa)), dims[0]*dims[1]*sizeof(double));
      for(i=0; i<par.N; ++i) {
        psol->x[i] = X[3*i];
        psol->y[i] = X[3*i+1];
        psol->z[i] = X[3*i+2];
      }
      free(X);
    } else if(strcmp(name, "Gamma") == 0 ) {
      memcpy((void *)pdyn->G, (void *)(mxGetPr(pa)), par.N*sizeof(double));
//       printf("Copied %d vortex strengths, of which the last one was %f\n", par.N, pdyn->G[par.N-1]);
    }
    mxDestroyArray(pa);
  }
  if (matClose(pmat) != 0) {
    printf("Error closing file %s\n",file);
    return(1);
  }
//   printf("Done reading array\n");
  return(0);
}

int readOrderingFromMatLab(struct Dynamics *pdyn, struct Parameters par, char *file)
{  
  MATFile *pmat;
  const char *name;
  int   ndim;
  int   i;
  mxArray *pa;
  const int *dims;
  
  dims = (int *) malloc(2*(sizeof(int)));
  
  printf("Reading file %s...\n\n", file);

  pmat = matOpen(file, "r");
  if (pmat == NULL) {
    printf("Error reopening file %s\n", file);
    return(1);
  }

  /* Get headers of all variables */
//   printf("\nExamining the header for each variable:\n");
  while( (pa = matGetNextVariable(pmat, &name)) != NULL )
  {
    ndim = mxGetNumberOfDimensions(pa);
    dims = mxGetDimensions(pa);
    if( strcmp(name, "A") == 0 ) {
      memcpy((void *)pdyn->A, (void *)(mxGetData(pa)), dims[0]*dims[1]*sizeof(unsigned char));
    } else if(strcmp(name, "B") == 0 ) {
      memcpy((void *)pdyn->B, (void *)(mxGetData(pa)), dims[0]*dims[1]*sizeof(unsigned char));
    }
    mxDestroyArray(pa);
  }

  if (matClose(pmat) != 0) {
    printf("Error closing file %s\n",file);
    return(1);
  }
//   printf("Done reading array\n");
  return(0);
}

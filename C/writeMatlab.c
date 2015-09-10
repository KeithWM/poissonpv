#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mat.h>
#if 0
#include <fftw3.h>
#endif

void writeArrayToMatlab(MATFile *pmatfile, double *store, char *name, int m, int n)
{
  mxArray *matlabArray;
  
  matlabArray = (mxArray *) mxCreateDoubleMatrix(m,n,mxREAL); /* creates a mxn MATLAB array of reals */
  memcpy((void *)(mxGetPr(matlabArray)),   (void *)store, m*n*sizeof(double)); /* copies the memory (m*n) at Q to tseries */
  matPutVariable(pmatfile,name,matlabArray); /* puts the array matlabx into the .mat file as "dt" */
  mxDestroyArray(matlabArray); /* clears the array (for memory I presume) */  
//   printf("Written %s to Matlabfile, it was %d by %d\n", name, m, n);
}

#if 0
void writeComplexArrayToMatlab(MATFile *pmatfile, fftw_complex *store, char *name, int m, int n)
{
  int i;
  
  mxArray *matlabArray;
  double *storeR, *storeI;
  
  storeR = (double *) malloc(m*n*sizeof(double));
  storeI = (double *) malloc(m*n*sizeof(double));
  
  for(i=0;i<m*n;++i) {
    storeR[i] = store[i][0];
    storeI[i] = store[i][1];
  }
  
  matlabArray = (mxArray *) mxCreateDoubleMatrix(m,n,mxCOMPLEX); /* creates a mxn MATLAB array of reals */
  memcpy((void *)(mxGetPr(matlabArray)),   (void *)storeR, m*n*sizeof(double)); /* copies the memory (m*n) at Q to tseries */
  memcpy((void *)(mxGetPi(matlabArray)),   (void *)storeI, m*n*sizeof(double)); /* copies the memory (m*n) at Q to tseries */
  matPutVariable(pmatfile,name,matlabArray); /* puts the array matlabx into the .mat file as "dt" */
  mxDestroyArray(matlabArray); /* clears the array (for memory I presume) */
  free(storeR);
  free(storeI);
  printf("Written %s to Matlabfile, it was %d by %d\n", name, m, n);
}
#endif 

void writeIntArrayToMatlab(MATFile *pmatfile, unsigned char *store, char *name, int m, int n)
{
  mxArray *matlabArray;
  
  matlabArray = (mxArray *) mxCreateNumericMatrix(m,n,mxINT32_CLASS,mxREAL); /* creates a mxn MATLAB array of reals */
  memcpy((void *)(mxGetData(matlabArray)),   (void *)store, m*n*sizeof(unsigned char)); /* copies the memory (m*n) at Q to tseries */
  matPutVariable(pmatfile,name,matlabArray); /* puts the array matlabx into the .mat file as "dt" */
  mxDestroyArray(matlabArray); /* clears the array (for memory I presume) */  
//   printf("Written %s to Matlabfile, it was %d by %d\n", name, m, n);
}

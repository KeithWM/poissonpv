#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mat.h>
#if 0
#include <fftw3.h>
#endif

void writeArrayToMatlab(MATFile *pmatfile, double *store, char *name, int m, int n);
//void writeComplexArrayToMatlab(MATFile *pmatfile, fftw_complex *store, char *name, int m, int n);
void writeIntArrayToMatlab(MATFile *pmatfile, unsigned char *store, char *name, int m, int n);

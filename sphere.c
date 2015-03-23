#define _REENTRANT
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <string.h>
#include "mat.h"
#include "normal.h"
#include <dispatch/dispatch.h>
#include <sys/time.h>
#include <sys/stat.h>
#include "readMatlab.h"
#include "writeMatlab.h"
#include "parallel.h"
#include "normal.h"

dispatch_semaphore_t  sema[THREADS*THREADS];
dispatch_semaphore_t gsema[2*THREADS*(4*THREADS-2)];

void initializeSol(struct Solution *psol, struct Parameters par) {
  psol->t = 0;
  psol->x = (double *) malloc(par.N*sizeof(double));
  psol->y = (double *) malloc(par.N*sizeof(double));
  psol->z = (double *) malloc(par.N*sizeof(double));
  
  psol->H = (double *) malloc(sizeof(double));
  psol->J = (double *) malloc(4*sizeof(double));
  psol->Hc = (double *) malloc(THREADS*sizeof(double));
  psol->Jc = (double *) malloc(THREADS*4*sizeof(double));
}

void initializeDyn(struct Dynamics *pdyn, struct Parameters par) {
  pdyn->G = (double *) malloc(par.N*sizeof(double));
  
  pdyn->A = (unsigned char *) malloc(THREADS*2*pdyn->R*sizeof(long int));
  pdyn->B = (unsigned char *) malloc(THREADS*2*pdyn->R*sizeof(long int));
}

void initializeSto(struct Storage *psto, struct Parameters par) {
  psto->t = (double *) malloc((par.Nouter+1)*sizeof(double));
  
  psto->x = (double *) malloc((par.Nouter+1)*par.N*sizeof(double));
  psto->y = (double *) malloc((par.Nouter+1)*par.N*sizeof(double));
  psto->z = (double *) malloc((par.Nouter+1)*par.N*sizeof(double));
  
  psto->H = (double *) malloc((par.Nouter+1)*sizeof(double));
  psto->J = (double *) malloc((par.Nouter+1)*4*sizeof(double));
}

int main(int argc, char *argv[]) {
  int i,j, t;
  struct Parameters par;
  struct Solution   sol;
  struct Dynamics   dyn;
  struct Storage    sto;
  struct timeval tv1, tv2;
  double *times;
  char orderfile[999];
  char icfile[999];
  char outputdir[999];
  char outputfile[999];
  void *Pass[4];
  MATFile *pmatfile;
  dispatch_queue_t queue;
  
  Pass[0] = &par;
  Pass[1] = &sol;
  Pass[2] = &dyn;
  Pass[3] = &sto;
  
  /* specify some defaults */
  par.N = 16*THREADS;
  par.T = 1;
  par.compo = "LT";
  par.order = "naive";
  par.energy= "neutral";
  par.Tend = 1;
  par.dt  = .1;
  par.dto = .1;
  par.synch = 1;
  strcpy(outputdir, ".");
  
  for (i=1; i<argc; ++i) {
    if(!strcmp(argv[i], "--help") || !strcmp(argv[i], "-h")) {
      printf(
"usage: %s [options] \n"
"where options are:\n"
"-N      specify number of vortices\n"
"-T      specify number of times to repeat experiment\n"
"-compo  specify type of composition LT or S (Lie-Trotter or Strang)\n"
"-order  specify type of ordering naive, smart or fractal \n"
"-dt     specify time step \n"
"-dto    specify output time step \n"
"-Tend   specify end time \n"
"-synch  specify synchronous simulation (1) or no (0)\n"
"-batch  specify a name for the batch\n"
"-energy specify energy level\n"
,argv[0]);
      return 1;
    }
    if (!strcmp(argv[i], "-N"))
      par.N = atoi(argv[++i]);
    if (!strcmp(argv[i], "-T"))
      par.T = atoi(argv[++i]);
    if (!strcmp(argv[i], "-compo")) {
      par.compo = argv[++i];
      if (strcmp(par.compo, "LT") && strcmp(par.compo, "S") && strcmp(par.compo, "Y4") && strcmp(par.compo, "M4") && strcmp(par.compo, "Y6") && strcmp(par.compo, "M6") ) {
        printf("Invalid composition specified for %s: %s\n", argv[i-1], argv[i]);
        return 1;
      }
    }
    if (!strcmp(argv[i], "-order")) {
      par.order = argv[++i];
      if (strcmp(par.order, "naive") && strcmp(par.order, "smart") && strcmp(par.order, "fractal") && strcmp(par.order, "fractal_new")) {
        printf("Invalid order specified for %s: %s\n", argv[i-1], argv[i]);
        return 1;
      }
    }
    if (!strcmp(argv[i], "-energy")) {
      par.energy = argv[++i];
    }
    if (!strcmp(argv[i], "-dt"))
      par.dt   = atof(argv[++i]);
    if (!strcmp(argv[i], "-dto"))
      par.dto  = atof(argv[++i]);
    if (!strcmp(argv[i], "-Tend"))
      par.Tend = atof(argv[++i]);
    if (!strcmp(argv[i], "-synch"))
      par.synch = atoi(argv[++i]);
    if (!strcmp(argv[i], "-batch"))
      strcpy(outputdir, argv[++i]);
  }
  times = (double *) malloc(par.T*sizeof(double));

  par.L = (par.N-1)/(2*THREADS)+1; // number of vortices per group
  par.Ng= 2*THREADS*par.L; // number of vortices, including "ghost" vortices
  par.Nouter = (int) round(par.Tend/par.dto);
  par.dto = (par.Nouter == 0) ? 1. : par.Tend/par.Nouter;
  par.Ninner = (int) round(par.dto/par.dt);
  par.Ninner =  par.Ninner;
  par.dt = (par.Ninner == 0) ? 1. : par.dto/par.Ninner;
  if(strcmp(par.compo, "LT") == 0 || strcmp(par.compo, "S") == 0) { // see Hairer, Lubich and Wanner for more details
    par.dt_split = par.dt; // default time step of a whole composition, will be overwritten for higher-order compositions, initialized here for LT & S for speed
  } else {
    if (strcmp(par.compo, "Y4") == 0) { // composition of THREE Strang Splittings for order 4 (Creutz & Gocksch 1989, Forest 1989, Yoshida 1990)
      par.Nsplit = 3;
      par.dts_split = (double *) malloc(par.Nsplit * sizeof(double));
      par.dts_split[0] =  1.351207191959658;
    } else if (strcmp(par.compo, "M4") == 0) { // composition of FIVE Strang Splittings for order 4 (McLachlan 1995)
      par.Nsplit = 5;
      par.dts_split = (double *) malloc(par.Nsplit * sizeof(double));
      par.dts_split[0] =  0.28;
      par.dts_split[1] =  0.62546642846767004501;
    } else if (strcmp(par.compo, "Y6") == 0) { // composition of SEVEN Strang Splittings for order 6 (Yoshida 1990)
      par.Nsplit = 7;
      par.dts_split = (double *) malloc(par.Nsplit * sizeof(double));
      par.dts_split[0] =  0.78451361047755726382;
      par.dts_split[1] =  0.23557321335935813368;
      par.dts_split[2] = -1.17767998417887100695;
    } else if (strcmp(par.compo, "M6") == 0) { // composition of NINE Strang Splittings for order 6 (McLachlan 1995)
      par.Nsplit = 9;
      par.dts_split = (double *) malloc(par.Nsplit * sizeof(double));
      par.dts_split[0] =  0.39216144400731413928;
      par.dts_split[1] =  0.33259913678935943860;
      par.dts_split[2] = -0.70624617255763935981;
      par.dts_split[3] =  0.08221359629355080023;
    }
    par.dts_split[(par.Nsplit-1)/2] = 1.;
    for (i=0; i<(par.Nsplit-1)/2; ++i) {
      par.dts_split[(par.Nsplit-1)/2]-= 2*(par.dts_split[par.Nsplit-1-i] = par.dts_split[i]);
    }
//    for (i=0; i<par.Nsplit; ++i) {
//      printf("w[%d]: %f\n", i, par.dts_split[i]);
//    }
  }
  printf("Using %d inner time steps of %f each, resulting in %d outer times steps of %f each.\n", par.Ninner, par.dt, par.Nouter, par.dto);
  
  dyn.R = 2*THREADS-1;
  dyn.S = !strcmp(par.compo, "LT") ? 2*THREADS-1 : 4*THREADS-2;
  
  sprintf(orderfile, "../orders/order_%s_C%d.mat", par.order, THREADS);
  sprintf(icfile, "../ic/ics_%s_M8N%d.mat", par.energy, par.N);
  mkdir(outputdir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  sprintf(outputdir, "%s/%s", outputdir, par.energy);
  mkdir(outputdir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  sprintf(outputdir, "%s/Tend%.2e", outputdir, par.Tend);
  mkdir(outputdir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  sprintf(outputdir, "%s/C%d_compo%s_synch%d_N%d_order%s_dt%.2e", outputdir, THREADS, par.compo, par.synch, par.N, par.order, par.dt);
  mkdir(outputdir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  
  initializeSol(&sol, par);
  initializeDyn(&dyn, par);
  initializeSto(&sto, par);
  readOrderingFromMatLab(&dyn, orderfile);

/*
  for(i=0; i<THREADS; ++i) {
    for(j=0; j<2*dyn.R; ++j) {
      printf("(%d,%d)\t", dyn.A[2*dyn.R*i + j], dyn.B[2*dyn.R*i + j]);
    }
    printf("\n");
  }
*/
//   printf("Using %d steps\n\n", dyn.S);
  
  for(t=0; t<par.T; ++t) { 
//    printf("Simulating a system of %d vortices on %d cores, with %d vortices per group, resulting in %d \"ghost\" vortices\n", par.N, THREADS, par.L, par.Ng-par.N);

    readInitCondFromMatLab(&sol, &dyn, par, icfile);
//    for(i=0; i<par.N; ++i) {
//      printf("G[%d]: %f, x[%d]: (%f,%f,%f)\n", i, dyn.G[i], i, sol.x[i], sol.y[i], sol.z[i]);
//    }
     
    queue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0);
    
    for(i=0; i<THREADS*THREADS; ++i) {
      sema[i] = dispatch_semaphore_create(0);
    }
    for(i=0; i<2*THREADS*(4*THREADS-2); ++i) {
      gsema[i] = dispatch_semaphore_create(0);
    }

    gettimeofday(&tv1, NULL);
    dispatch_apply_f(THREADS, queue, Pass, (void *) parallel);
    gettimeofday(&tv2, NULL);
    
    times[t] = (double) (tv2.tv_usec - tv1.tv_usec)/1000000 + (double) (tv2.tv_sec - tv1.tv_sec);
    printf ("Total time = %f seconds\n", times[t]);
    
    sprintf(outputfile, "%s/expid%0ld%0d.mat", outputdir, tv2.tv_sec, tv2.tv_usec);
//    sprintf(outputfile, "test.mat");
    if((pmatfile = matOpen(outputfile,"w"))) {
      writeArrayToMatlab(   pmatfile, dyn.G,  "G",    par.N,   1);
      writeIntArrayToMatlab(pmatfile, dyn.A,  "A",    THREADS, 2*(2*THREADS-1));
      writeIntArrayToMatlab(pmatfile, dyn.B,  "B",    THREADS, 2*(2*THREADS-1));
      
      writeArrayToMatlab(pmatfile, times+t,"time", 1,       1);
      
      writeArrayToMatlab(pmatfile, sto.t,"t",      1, par.Nouter+1);
      if(par.N < 1000) {
        writeArrayToMatlab(pmatfile, sto.x,"x",      par.N, par.Nouter+1);
        writeArrayToMatlab(pmatfile, sto.y,"y",      par.N, par.Nouter+1);
        writeArrayToMatlab(pmatfile, sto.z,"z",      par.N, par.Nouter+1);
      }

      writeArrayToMatlab(pmatfile, sto.H,"H",      1, par.Nouter+1);
      writeArrayToMatlab(pmatfile, sto.J,"J",      4, par.Nouter+1);
      
      printf("MatLab file stored at %s\n", outputfile);
    } else {
      printf("failed to initalize pmatfile\n");
    }
  }
}

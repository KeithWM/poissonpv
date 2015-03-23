#ifndef PARALLEL_H
#define PARALLEL_H

struct Parameters {
  int N; // number of vortices
  int L; // number of vortices per group
  int Ng;// number of vortices, including "ghosts"
  int P; // THREADS (not used)
  int T; // number of times to repeat the experiment
  int Nouter; // number of steps in outer loop
  int Ninner; // number of steps in inner loop
  int synch; // synchronour solution or not

  char *compo;
  char *order;
  char *energy;

  double dt;
  double dt_split;
  double *dts_split;
  int    Nsplit;
  double dto;
  double Tend;
};

struct Solution {
  double t;
  
  double *x;
  double *y;
  double *z;
  
  double *H;
  double *J;
  
  double *Hc;
  double *Jc;
};

struct Dynamics {
  unsigned char *A;
  unsigned char *B;
  //long int *A;  
  //long int *B;  
  int R; // number of rounds (in Lie-Trotter terms, 2R in Strang)
  int S; // number of rounds actually in this scheme (in Lie-Trotter terms, 2R in Strang)
  
  double *G;
};

struct Storage {
  double *t;
  
  double *x;
  double *y;
  double *z;
  
  double *H;
  double *J;
};

void parallel(void **Pass, size_t thread);
int pairbetween(int *i, int *j, int counter, int L, int dir);
int pairwithin(int *i, int *j, int counter, int L, int dir);
void waitfor(int k, int l, size_t thread, int s);
void donewith(int k, int l, size_t thread, int s);
void barrier(size_t thread);

// spherestep.c
void spherestep(struct Solution *psol, struct Parameters parc, struct Dynamics dync, size_t thread);

// conserved.c
// double Hpair(int i, int j, struct Solution *psol, struct Parameters parc, struct Dynamics dync);
// void incHbetween(double *Hct, unsigned char k, unsigned char l, struct Solution *psol, struct Parameters parc, struct Dynamics dync);
// void incHwithin(double *Hct, unsigned char k, unsigned char l, struct Solution *psol, struct Parameters parc, struct Dynamics dync);
// void incJ(double *Jct, int k, struct Solution *psol, struct Parameters parc, struct Dynamics dync);
void conserved(struct Solution *psol, struct Parameters parc, struct Dynamics dync, size_t thread);
// store.c
void store(struct Solution *psol, struct Storage *psto, struct Parameters parc, size_t thread, int no);
#endif

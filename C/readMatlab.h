#ifndef READMATLAB_H
#define READMATLAB_H

#include "parallel.h"

int readInitCondFromMatLab(struct Solution *psol, struct Dynamics *pdyn, struct Parameters par, char *file);
int readOrderingFromMatLab(struct Dynamics *pdyn, char *file);
#endif


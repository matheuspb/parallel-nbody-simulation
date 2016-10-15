#ifndef PARALLEL_H
#define PARALLEL_H

#include <pthread.h>

void ComputeForces();
void* ComputeForcesThread(void*);
void ComputeNewPos();
void* ComputeNewPosThread(void*);
void InitParticles();
double Random(void);

typedef struct {
	double x, y, z;
	double mass;
} Particle;

typedef struct {
	double xold, yold, zold;
	double fx, fy, fz;
} ParticleV;

#endif


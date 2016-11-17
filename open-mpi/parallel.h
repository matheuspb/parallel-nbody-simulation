#ifndef PARALLEL_H
#define PARALLEL_H

void ComputeForces();
void ComputeForcesThread();
void ComputeNewPos();
void ComputeNewPosThread();
void InitParticles();

typedef struct {
	double x, y, z;
	double mass;
} Particle;

typedef struct {
	double xold, yold, zold;
	double fx, fy, fz;
} ParticleV;

#endif


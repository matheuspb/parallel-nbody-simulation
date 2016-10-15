#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "parallel.h"

#define N_THREADS 4

double dt, dt_old;
int npart; // number of particles
Particle  * particles;   /* Particles */
ParticleV * pv;          /* Particle velocity */
double max_f, sim_t, chunk_size;
pthread_t compute_forces[N_THREADS];
pthread_t compute_pos[N_THREADS];
pthread_mutex_t max_lock;

int main(int argc, char **argv)
{
	int cnt; /* number of times in loop */

	if(argc != 3){
		printf("Wrong number of parameters.\n");
		printf("Usage: nbody num_bodies timesteps\n");
		exit(1);
	}

	npart = atoi(argv[1]);
	cnt = atoi(argv[2]);
	dt = 0.001;
	dt_old = 0.001;
	chunk_size = (double)npart/(double)N_THREADS;

	/* Allocate memory for particles */
	particles = (Particle *) malloc(sizeof(Particle)*npart);
	pv = (ParticleV *) malloc(sizeof(ParticleV)*npart);

	/* Generate the initial values */
	InitParticles();
	sim_t = 0.0;

	pthread_mutex_init(&max_lock, NULL);

	while (cnt--) {
		/* Compute forces (2D only) */
		ComputeForces();
		/* Once we have the forces, we compute the changes in position */
		ComputeNewPos();
	}

	pthread_mutex_destroy(&max_lock);

	//for (int i = 0; i < npart; i++)
	//	fprintf(stdout,"%.5lf %.5lf\n", particles[i].x, particles[i].y);

	free(particles);
	free(pv);

	return 0;
}

void InitParticles()
{
	int i;
	for (i=0; i<npart; i++) {
		particles[i].x	  = Random();
		particles[i].y	  = Random();
		particles[i].z	  = Random();
		particles[i].mass = 1.0;
		pv[i].xold	  = particles[i].x;
		pv[i].yold	  = particles[i].y;
		pv[i].zold	  = particles[i].z;
		pv[i].fx	  = 0;
		pv[i].fy	  = 0;
		pv[i].fz	  = 0;
	}
}

void ComputeForces()
{
	max_f = 0.0;
	int indexes[N_THREADS];
	for (int i = 0; i < N_THREADS; ++i) {
		indexes[i] = i;
		pthread_create(&compute_forces[i], NULL, ComputeForcesThread,
				(void*) &indexes[i]);
	}
	for (int i = 0; i < N_THREADS; ++i) {
		pthread_join(compute_forces[i], NULL);
	}
}

void* ComputeForcesThread(void* a) {
	int nt = *((int*)a);
	int begin = nt*chunk_size;
	int end = (nt+1)*chunk_size;
	for (int i=begin; i<end; i++) {
		int j;
		double xi, yi, rx, ry, mj, r, fx, fy, rmin;
		rmin = 100.0;
		xi   = particles[i].x;
		yi   = particles[i].y;
		fx   = 0.0;
		fy   = 0.0;
		for (j=0; j<npart; j++) {
			rx = xi - particles[j].x;
			ry = yi - particles[j].y;
			mj = particles[j].mass;
			r  = rx * rx + ry * ry;
			/* ignore overlap and same particle */
			if (r == 0.0) continue;
			if (r < rmin) rmin = r;
			r  = r * sqrt(r);
			fx -= mj * rx / r;
			fy -= mj * ry / r;
		}
		pv[i].fx += fx;
		pv[i].fy += fy;
		fx = sqrt(fx*fx + fy*fy)/rmin;
		pthread_mutex_lock(&max_lock);
		if (fx > max_f) max_f = fx;
		pthread_mutex_unlock(&max_lock);
	}
	return NULL;
}

double a0, a1, a2;

void ComputeNewPos()
{
	double dt_new;
	a0	 = 2.0 / (dt * (dt + dt_old));
	a2	 = 2.0 / (dt_old * (dt + dt_old));
	a1	 = -(a0 + a2);
	int indexes[N_THREADS];
	for (int i = 0; i < N_THREADS; ++i) {
		indexes[i] = i;
		pthread_create(&compute_pos[i], NULL, ComputeNewPosThread,
				(void*) &indexes[i]);
	}
	for (int i = 0; i < N_THREADS; ++i) {
		pthread_join(compute_pos[i], NULL);
	}
	dt_new = 1.0/sqrt(max_f);
	/* Set a minimum: */
	if (dt_new < 1.0e-6) dt_new = 1.0e-6;
	/* Modify time step */
	if (dt_new < dt) {
		dt_old = dt;
		dt     = dt_new;
	}
	else if (dt_new > 4.0 * dt) {
		dt_old = dt;
		dt    *= 2.0;
	}
	sim_t += dt_old;
}

void* ComputeNewPosThread(void* a) {
	int nt = *((int*)a);
	int begin = nt*chunk_size;
	int end = (nt+1)*chunk_size;
	for (int i = begin; i < end; ++i) {
		double xi, yi;
		xi	           = particles[i].x;
		yi	           = particles[i].y;
		particles[i].x = (pv[i].fx - a1 * xi - a2 * pv[i].xold) / a0;
		particles[i].y = (pv[i].fy - a1 * yi - a2 * pv[i].yold) / a0;
		pv[i].xold     = xi;
		pv[i].yold     = yi;
		pv[i].fx       = 0;
		pv[i].fy       = 0;
	}
	return NULL;
}

/*
pRNG based on http://www.cs.wm.edu/~va/software/park/park.html
*/
#define MODULUS    2147483647
#define MULTIPLIER 48271
#define DEFAULT    123456789

static long seed = DEFAULT;

double Random(void)
/*
----------------------------------------------------------------
Random returns a pseudo-random real number uniformly distributed
between 0.0 and 1.0.
----------------------------------------------------------------
*/
{
	const long Q = MODULUS / MULTIPLIER;
	const long R = MODULUS % MULTIPLIER;
	long t;

	t = MULTIPLIER * (seed % Q) - R * (seed / Q);
	if (t > 0)
		seed = t;
	else
		seed = t + MODULUS;
	return ((double) seed / MODULUS);
}

/*
End of the pRNG algorithm
*/


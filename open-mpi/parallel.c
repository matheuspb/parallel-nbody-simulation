#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "parallel.h"
#include "random.h"

int npart; // number of particles
double max_f, sim_t, dt, dt_old;
double chunk_size; // number of particles processed by each process

double a0, a1, a2; // used in compute new position functions

int size, rank;

Particle  * particles;   /* Particles */
ParticleV * pv;          /* Particle velocity */

int main(int argc, char **argv) {
	//MPI_Status st;
	MPI_Init(&argc,	&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == 0) {
		printf("-np %d\n", size);
		if(argc != 4){
			printf("Wrong number of parameters.\n");
			printf("Usage: %s num_bodies timesteps print_results(0, !=0)\n"
					,argv[0]);
			exit(1);
		}

		npart = atoi(argv[1]);
		int cnt = atoi(argv[2]); /* number of times in loop */

		dt = 0.001;
		dt_old = 0.001;

		chunk_size = (double)npart;// /(double)(size - 1);

		/* Allocate memory for particles */
		particles = (Particle *) malloc(sizeof(Particle)*npart);
		pv = (ParticleV *) malloc(sizeof(ParticleV)*npart);

		sim_t = 0.0;
		InitParticles();

		while (cnt--) {
			ComputeForces();
			ComputeNewPos();
		}

		if (atoi(argv[3])) {
			for (int i = 0; i < npart; i++)
				fprintf(stdout,"%.5lf %.5lf\n", particles[i].x, particles[i].y);
		}

		free(particles);
		free(pv);
	}

	MPI_Finalize();
	return 0;
}

void InitParticles() {
	int i;
	for (i = 0; i < npart; i++) {
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

void ComputeForces() {
	max_f = 0.0;
	ComputeForcesThread();
}

void ComputeForcesThread() {
	int nt = rank;
	int begin = nt*chunk_size;
	int end = (nt+1)*chunk_size;

	for (int i = begin; i < end; i++) {
		double xi, yi, rx, ry, mj, r, fx, fy, rmin;
		rmin = 100.0;

		xi   = particles[i].x;
		yi   = particles[i].y;

		fx   = 0.0;
		fy   = 0.0;

		for (int j = 0; j < npart; j++) {
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

		if (fx > max_f) max_f = fx;
	}
}


void ComputeNewPos() {
	double dt_new;
	a0	 = 2.0 / (dt * (dt + dt_old));
	a2	 = 2.0 / (dt_old * (dt + dt_old));
	a1	 = -(a0 + a2);

	ComputeNewPosThread();

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

void ComputeNewPosThread() {
	int nt = rank;
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
}

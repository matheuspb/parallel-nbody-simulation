.PHONY: parallel open-mpi

all: nbody parallel open-mpi

nbody: nbody.c
	gcc -o nbody.out nbody.c -std=c11 -lm

parallel: parallel/parallel.c parallel/random.c
	gcc -o parallel.out -std=c11 -lpthread -lm -Wall parallel/*.c

open-mpi: open-mpi/parallel.c open-mpi/random.c
	mpicc -Wall -lm -o open-mpi.out open-mpi/*.c

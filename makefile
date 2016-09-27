all: nbody.c parallel.c
	gcc -o nbody.out nbody.c -std=c11 -lm
	gcc -o parallel.out parallel.c -std=c11 -lpthread -lm -Wall

nbody: nbody.c
	gcc -o nbody.out nbody.c -std=c11 -lm

parallel: parallel.c
	gcc -o parallel.out parallel.c -std=c11 -lpthread -lm -Wall


all: nbody.c parallel.c
	gcc -o nbody.out nbody.c -std=c11 -lm
	gcc -o parallel.out parallel.c -std=c11 -lpthread -lm -Wall


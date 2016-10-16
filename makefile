all: nbody parallel

nbody: nbody.c
	gcc -o nbody.out nbody.c -std=c11 -lm

parallel: parallel.c random.c
	gcc -o parallel.out -std=c11 -lpthread -lm -Wall parallel.c random.c


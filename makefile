.PHONY: parallel

CFLAGS = -std=c11 -lm

default: nbody parallel
	@echo "Succesfully compiled with $(CC)"

nbody: nbody.c
	@$(CC) $(CFLAGS) nbody.c -o nbody.out

parallel: parallel/parallel.c parallel/random.c
	@$(CC) $(CFLAGS) -lpthread -Wall parallel/*.c -o parallel.out

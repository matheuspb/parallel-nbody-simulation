.PHONY: end clean

CFLAGS = -std=c11 -pedantic -march=native -O3
LDFLAGS = -lm

SRC = $(wildcard src/*.c)
OBJ = $(SRC:.c=.o)

all: nbody parallel

nbody: $(nbody.c)

parallel: src/parallel
	cp src/parallel .

src/parallel: CFLAGS += -Wall -Wextra -Werror
src/parallel: LDFLAGS += -lpthread
src/parallel: $(OBJ)

clean:
	@rm -f $(OBJ)
	@rm -f nbody
	@rm -f parallel
	@rm -f src/parallel

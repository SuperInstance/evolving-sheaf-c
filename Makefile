CC = gcc
CFLAGS = -O2 -Wall -Wextra -I include
LDFLAGS = -lm

SRC = src/evolving_sheaf.c
TEST_SRC = tests/test_evolving_sheaf.c

.PHONY: all test clean

all: test_runner

test_runner: $(TEST_SRC) $(SRC) include/evolving_sheaf.h
	$(CC) $(CFLAGS) -o $@ $(TEST_SRC) $(SRC) $(LDFLAGS)

test: test_runner
	./test_runner

clean:
	rm -f test_runner

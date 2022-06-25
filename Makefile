# Makefile for the genomikon

LIB = -lm
CC = gcc
CFLAGS = -O2 -Wall -Werror

OBJECTS = \
	align.o\
	feature.o\
	model.o\
	sequence.o\
	toolbox.o\

ARC = libgenomikon.a

APP = testing
SRC = testing.c
OBJ = testing.o

# Targets

default:
	make $(ARC)
	make $(APP)
	cd dusty && make
	cd hmmstar && make
	cd isoformer && make
	cd presti && make
	cd smithy && make
	cd wordy && make

$(ARC): $(OBJECTS)
	ar rvs $(ARC) $(OBJECTS)

$(APP): $(OBJ) $(OBJECTS)
	$(CC) -o $(APP) $(OBJ) $(OBJECTS) $(LIB)

clean:
	rm -f *.o $(APP) $(ARC)
	cd dusty && make clean
	cd hmmstar && make clean
	cd isoformer && make clean
	cd presti && make clean
	cd smithy && make clean
	cd wordy && make clean

test: $(ARC)
	make
	cd dusty && make test
	cd hmmstar && make test
	cd isoformer && make test
	cd presti && make test
	cd smithy && make test
	cd wordy && make test

depend: $(OBJECTS:.o=.c)
	gcc -MM $^ > $@

# Inference Rules

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

# Dependancies

include depend

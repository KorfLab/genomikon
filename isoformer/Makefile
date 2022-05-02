# Makefile for isoformer

LIB = -L.. -lgenomikon -lm
INC = -I..
CC = gcc
CFLAGS = -O2 -Wall -Werror

OBJECTS = isoform.o

APP1 = isoformer
OBJ1 = isoformer.o

APP2 = isocounter
OBJ2 = isocounter.o

APP3 = isorandom
OBJ3 = isorandom.o

default:
	make $(APP1)
	make $(APP2)
	make $(APP3)

$(APP1): $(OBJ1) $(OBJECTS)
	$(CC) -o $(APP1) $(OBJ1) $(OBJECTS) $(LIB)

$(APP2): $(OBJ2) $(OBJECTS)
	$(CC) -o $(APP2) $(OBJ2) $(OBJECTS) $(LIB)

$(APP3): $(OBJ3) $(OBJECTS)
	$(CC) -o $(APP3) $(OBJ3) $(OBJECTS) $(LIB)

test: $(APP1)
	perl test.pl

clean:
	rm -f *.o $(APP1) $(APP2) $(APP3)

%.o: %.c
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<

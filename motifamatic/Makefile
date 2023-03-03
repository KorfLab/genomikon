# Makefile for motifamatic and associated programs and libraries

LIB = -L.. -lgenomikon -lm
INC = -I..
CC = gcc
CFLAGS = -O2 -Wall -Werror

OBJECTS = dmg.o pwm.o

APP1 = testing
OBJ1 = testing.o

APP2 = motifamatic
OBJ2 = motifamatic.o

default:
	make $(APP1)
	make $(APP2)

$(APP1): $(OBJ1) $(OBJECTS)
	$(CC) -o $(APP1) $(OBJ1) $(OBJECTS) $(LIB)

$(APP2): $(OBJ2) $(OBJECTS)
	$(CC) -o $(APP2) $(OBJ2) $(OBJECTS) $(LIB)

test: $(APP1)
	perl test.pl

clean:
	rm -f *.o $(APP1) $(APP2) $(APP3)

%.o: %.c
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<

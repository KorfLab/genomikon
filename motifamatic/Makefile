# Makefile for motifamatic

LIB = -L.. -lgenomikon -lm
INC = -I..
CC = gcc
CFLAGS = -O2 -Wall -Werror

OBJECTS =

APP1 = motifamatic
OBJ1 = motifamatic.o

default:
	make $(APP1)

$(APP1): $(OBJ1) $(OBJECTS)
	$(CC) -o $(APP1) $(OBJ1) $(OBJECTS) $(LIB)

test: $(APP1)

clean:
	rm -f *.o $(APP1) $(APP2) $(APP3)

%.o: %.c
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<

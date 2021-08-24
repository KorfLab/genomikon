# Makefile for genomikon demo programs

LIB = -lm -lgenomikon
LINC = -L..
INC = -I..
CC = gcc
CFLAGS = -O2 -Wall -Werror

APP1 = dusty
OBJ1 = dusty.o

APP2 = geney
OBJ2 = geney.o

APP3 = smithy
OBJ3 = smithy.o

APP4 = viterby
OBJ4 = viterby.o

default:
	make $(APP1)
	make $(APP2)
	make $(APP3)
	make $(APP4)

$(APP1): $(OBJ1)
	$(CC) -o $(APP1) $(OBJ1) $(LINC) $(LIB)

$(APP2): $(OBJ2)
	$(CC) -o $(APP2) $(OBJ2) $(LINC) $(LIB)

$(APP3): $(OBJ3)
	$(CC) -o $(APP3) $(OBJ3) $(LINC) $(LIB)

$(APP4): $(OBJ4)
	$(CC) -o $(APP4) $(OBJ4) $(LINC) $(LIB)

clean:
	rm -f *.o $(APP1) $(APP2) $(APP3) $(APP4)

%.o: %.c
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<

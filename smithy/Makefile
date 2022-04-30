# Makefile for genomikon demo programs

LIB = -L.. -lgenomikon -lm
INC = -I..
CC = gcc
CFLAGS = -O2 -Wall -Werror

APP = smithy
OBJ = smithy.o

default:
	make $(APP)

$(APP): $(OBJ)
	$(CC) -o $(APP) $(OBJ) $(LIB)

clean:
	rm -f *.o $(APP)

test: $(APP)
	echo "tests needed..."

%.o: %.c
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<

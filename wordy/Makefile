# Makefile for genomikon demo programs

LIB = -L.. -lgenomikon -lm
INC = -I..
CC = gcc
CFLAGS = -O2 -Wall -Werror

APP = wordy
OBJ = wordy.o

default:
	make $(APP)

$(APP): $(OBJ)
	$(CC) -o $(APP) $(OBJ) $(LIB)

clean:
	rm -f *.o $(APP)

test: $(APP)
	echo "test me please"

%.o: %.c
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<

# Makefile for genomikon demo programs

LIB = -L.. -lgenomikon -lm
INC = -I..
CC = gcc
#CFLAGS = -O2 -Wall -Werror

APP1 = dusty
OBJ1 = dusty.o

APP2 = geney
OBJ2 = geney.o

APP3 = smithy
OBJ3 = smithy.o

APP4 = viterby
OBJ4 = viterby.o

APP5 = wordy
OBJ5 = wordy.o

default:
	make $(APP1)
	make $(APP2)
	make $(APP3)
	make $(APP4)
	make $(APP5)

$(APP1): $(OBJ1)
	$(CC) -o $(APP1) $(OBJ1) $(LIB)

$(APP2): $(OBJ2)
	$(CC) -o $(APP2) $(OBJ2) $(LIB)

$(APP3): $(OBJ3)
	$(CC) -o $(APP3) $(OBJ3) $(LIB)

$(APP4): $(OBJ4)
	$(CC) -o $(APP4) $(OBJ4) $(LIB)

$(APP5): $(OBJ5)
	$(CC) -o $(APP5) $(OBJ5) $(LIB)

clean:
	rm -f *.o $(APP1) $(APP2) $(APP3) $(APP4) $(APP5)

test: $(APP1) $(APP2) $(APP3) $(APP4) $(APP5)
	./dusty 777.fa
	./geney 777.fa 777.gff -apwm acceptor.pwm -dpwm donor.pwm -emm exon.mm -imm intron.mm -elen exon.len -ilen intron.len
	./smithy WHATEVER WHICHEVER
	./wordy popular.txt -pop 100 -gen 20 -seed 1

%.o: %.c
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<

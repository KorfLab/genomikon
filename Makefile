# Makefile for the genomikon

LIB = -lm
CC = gcc
CFLAGS = -O2 -Wall -Werror

OBJECTS = \
	align.o\
	feature.o\
	hmm.o\
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

$(ARC): $(OBJECTS)
	ar rvs $(ARC) $(OBJECTS)

$(APP): $(OBJ) $(OBJECTS)
	$(CC) -o $(APP) $(OBJ) $(OBJECTS) $(LIB)

clean:
	rm -f *.o $(APP) $(ARC)
	cd demo && make clean

depend: $(OBJECTS:.o=.c)
	gcc -MM $^ > $@

test: $(APP) $(ARC)
	./$(APP) -vec -ivec -fvec -tvec -map -tmap -smat -sw \
		-pipe Makefile -fasta demo/777.fa -gff demo/777.gff\
		-pwm demo/donor.pwm -mm demo/exon.mm -len demo/intron.len

demo: $(ARC) $(APP)
	cd demo && make

# Inference Rules

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

# Dependancies

include depend

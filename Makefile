##############################
# Makefile for the genomikon #
##############################

LIB = -lm
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

###########
# Targets #
###########

default:
	make $(ARC)
	make $(APP)

$(ARC): $(OBJECTS)
	ar rvs $(ARC) $(OBJECTS)

$(APP): $(OBJ) $(OBJECTS)
	$(CC) -o $(APP) $(OBJ) $(OBJECTS) $(LIB)

clean:
	rm -f *.o $(APP) $(ARC)
	cd dusty && make clean
	cd smithy && make clean
	cd geney && make clean

depend: $(OBJECTS:.o=.c)
	gcc -MM $^ > $@

test: $(APP) $(ARC)
	./$(APP) -vec -ivec -fvec -tvec -map -tmap -smat -sw \
		-pipe Makefile -fasta data/777.fa -gff data/777.gff\
		-pwm data/donor.pwm -mm data/exon.mm -len data/intron.len

all: $(ARC) $(APP)
	cd dusty && make
	cd smithy && make
	cd geney && make

###################
# Inference Rules #
###################

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

################
# Dependancies #
################

include depend

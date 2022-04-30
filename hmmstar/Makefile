##############################
# Makefile for Prestigitizer #
##############################

LIB = -L.. -lgenomikon -lm
INC = -I..
CC = gcc
CFLAGS = -O2 -Wall -Werror

#OBJECTS = \

APP = hmmstar
SRC = hmmstar.c
OBJ = hmmstar.o

###########
# Targets #
###########

default:
	make gcc

$(APP): $(OBJ) $(OBJECTS)
	$(CC) -o $(APP) $(OBJ) $(OBJECTS) $(LIB)

test:
	echo "insert test here"

clean:
	rm -f *.o $(APP)

#################
# Architectures #
#################

gcc:
	make $(APP)


###################
# Inference Rules #
###################

%.o: %.c
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<




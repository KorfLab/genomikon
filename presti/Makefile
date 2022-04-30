##############################
# Makefile for Prestigitizer #
##############################

LIB = -L.. -lgenomikon -lm
INC = -I..
CC = gcc
CFLAGS = -O2 -Wall -Werror

#OBJECTS = \

APP = presti
SRC = presti.c
OBJ = presti.o

###########
# Targets #
###########

default:
	make gcc

$(APP): $(OBJ) $(OBJECTS)
	$(CC) -o $(APP) $(OBJ) $(OBJECTS) $(LIB)

test: $(APP)
	./$(APP) values.txt 1 5 20
	echo "proper test needed"

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




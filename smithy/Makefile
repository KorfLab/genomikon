# Makefile for smithy

LIB = -lm -lgenomikon
LINC = -L..
INC = -I..
CFLAGS = -O2 -Wall -Werror

APP = smithy
OBJ = smithy.o

default:
	make $(APP)

$(APP): $(OBJ) $(OBJECTS)
	$(CC) -o $(APP) $(OBJ) $(OBJECTS) $(LINC) $(LIB)

clean:
	rm -f *.o $(APP)

%.o: %.c
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<

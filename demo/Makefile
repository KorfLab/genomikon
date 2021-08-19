# Makefile for genomikon demo programs

LIB = -lm -lgenomikon
LINC = -L..
INC = -I..
CFLAGS = -O2 -Wall -Werror

APP1 = dusty
OBJ1 = dusty.o

APP2 = geney
OBJ2 = geney.o

APP3 = smithy
OBJ3 = smithy.o

default:
	make $(APP1)
	make $(APP2)
	make $(APP3)

$(APP1): $(OBJ1)
	$(CC) -o $(APP1) $(OBJ1) $(LINC) $(LIB)

$(APP2): $(OBJ2)
	$(CC) -o $(APP2) $(OBJ2) $(LINC) $(LIB)

$(APP3): $(OBJ3)
	$(CC) -o $(APP3) $(OBJ3) $(LINC) $(LIB)

clean:
	rm -f *.o $(APP1) $(APP2) $(APP3)

%.o: %.c
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<

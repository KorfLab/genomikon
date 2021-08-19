# Makefile for geney

LIB = -lm -lgenomikon
LINC = -L..
INC = -I..
CFLAGS = -O2 -Wall -Werror

APP = geney
OBJ = geney.o

default:
	make $(APP)

$(APP): $(OBJ) $(OBJECTS)
	$(CC) -o $(APP) $(OBJ) $(OBJECTS) $(LINC) $(LIB)

test: $(APP)
	./$(APP) 777.fa 777.gff \
		-dpwm donor.pwm -apwm acceptor.pwm \
		-emm exon.mm -imm intron.mm \
		-elen exon.len -ilen intron.len

clean:
	rm -f *.o $(APP)

%.o: %.c
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<

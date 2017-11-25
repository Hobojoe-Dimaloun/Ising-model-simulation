IDIR =../include
CC=gcc
CFLAGS=-Wall -g -fopenmp -I$(IDIR)

ODIR=obj
LDIR =../lib

LIBS=-lm -lSDL2 -lSDL2_image -lgsl -lgslcblas

_DEPS =
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = main.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

../bin/executable: $(OBJ)
	gcc -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~

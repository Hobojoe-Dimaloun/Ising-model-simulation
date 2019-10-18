IDIR =../include
CC=gcc
CFLAGS=-Wall -g -fopenmp -I$(IDIR) -IC:/GSL/include

ODIR=obj
LDIR =../lib

LIBS=-lm -lgsl -lgslcblas -lmsmpi -L$(LDIR) -LC:/GSL/lib

_DEPS = ising_functions.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = ising_functions.o main.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

../bin/executable: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~

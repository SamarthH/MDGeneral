IDIR =../include
CC=g++
CFLAGS=-I$(IDIR) -fopenmp -Wall

ODIR=obj
LDIR =../lib

_DEPS = algorithm_constants.h client.h constants.h initialize.h integrate.h interaction.h system.h thermo.h thermostat.h write.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = client.o initialize.o integrate.o interaction.o thermo.o thermostat.o write.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

mdgen_back: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(IDIR)/*~ 

IDIR =../include
CC=g++
CFLAGS=-I$(IDIR) -Wall 

ODIR=obj
LDIR =../lib

LIBS= -ltrng4 -fopenmp

_DEPS = algorithm_constants.h client.h constants.h initialize.h integrate.h interaction.h system.h thermo.h thermostat.h write.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = client.o initialize.o integrate.o interaction.o thermo.o thermostat.o write.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.cpp $(DEPS)
	$(CC) -L$(LDIR) -c -o $@ $< $(CFLAGS) $(LIBS)

mdgen_back: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(IDIR)/*~ 

nest:
	set OMP_DYNAMIC=true
	set OMP_NESTED=true
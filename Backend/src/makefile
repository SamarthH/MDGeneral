IDIR =../include
CC=g++
CFLAGS=-I$(IDIR) -Wall 
EXTRAFLAGS=

ODIR=obj
LDIR =../lib

LIBS= -ltrng4 -fopenmp

_DEPS = algorithm_constants.h client.h constants.h initialize.h integrate.h interaction.h quaternion.h system.h thermo.h thermostat.h universal_functions.h write.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = client.o initialize.o integrate.o interaction.o quaternion.o system.o thermo.o thermostat.o write.o universal_functions.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.cpp $(DEPS)
	$(CC) -L$(LDIR) -c -o $@ $< $(CFLAGS) $(EXTRAFLAGS) $(LIBS)

mdgen_back: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(EXTRAFLAGS) $(LIBS)

	mv mdgen_back ../

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(IDIR)/*~ 

nest:
	set OMP_DYNAMIC=true
	set OMP_NESTED=true
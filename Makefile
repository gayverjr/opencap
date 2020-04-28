CXX=g++-9
CXXFLAGS= -Wall -g -std=gnu++17 -fopenmp
LIBS=-larmadillo -lm -lnumgrid -lhdf5 
DEPS = System.h Shell.h overlap.h Atom.h BasisSet.h utils.h transforms.h CAP.h \
molden_transform.h read_qchem_fchk.h read_rassi_h5.h gto_ordering.h InputParser.h BasisSetParser.h
OBJ = System.o Shell.o BasisSet.o Atom.o overlap.o main.o utils.o transforms.o \
read_qchem_fchk.o read_rassi_h5.o gto_ordering.o InputParser.o BasisSetParser.o CAP.o


%.o: %.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

main: $(OBJ) 
	$(CXX) $(CXXFLAGS) -o main $(OBJ) $(LIBS)

.PHONY: clean

clean:
	rm *.o
	

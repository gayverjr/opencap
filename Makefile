CXX=g++-9
CXXFLAGS= -Wall -g -std=gnu++17 -fopenmp
LIBS=-larmadillo -lm -lnumgrid -lhdf5 
DEPS = System.h Shell.h overlap.h Atom.h BasisSet.h utils.h transforms.h numerical.h CAP.h \
molden_transform.h read_qchem_fchk.h readMolcasHDF5.h molcas_transform.h InputParser.h BasisSetParser.h
OBJ = System.o Shell.o BasisSet.o Atom.o overlap.o main.o utils.o transforms.o numerical.o molden_transform.o \
read_qchem_fchk.o readMolcasHDF5.o molcas_transform.o InputParser.o BasisSetParser.o CAP.o


%.o: %.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

main: $(OBJ) 
	$(CXX) $(CXXFLAGS) -o main $(OBJ) $(LIBS)

.PHONY: clean

clean:
	rm *.o
	

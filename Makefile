CXX=g++-9
CXXFLAGS= -Wall -g -std=gnu++17 -fopenmp
LIBS=-larmadillo -lm -lnumgrid -lhdf5 
DEPS = Shell.h overlap.h Atom.h BasisSet.h utils.h transforms.h System.h numerical.h CAP.h \
molden_transform.h read_qchem_fchk.h readMolcasHDF5.h molcas_transform.h
OBJ = Shell.o BasisSet.o Atom.o overlap.o libCAP.o utils.o transforms.o System.o numerical.o molden_transform.o \
read_qchem_fchk.o readMolcasHDF5.o molcas_transform.o


%.o: %.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

main: $(OBJ) 
	$(CXX) $(CXXFLAGS) -o main $(OBJ) $(LIBS)

.PHONY: clean

clean:
	rm *.o
	

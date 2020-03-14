CXX=g++-9
CXXFLAGS= -Wall -g -std=gnu++17 -fopenmp
LIBS=-larmadillo -lm -lnumgrid
DEPS = Shell.h overlap.h Atom.h BasisSet.h utils.h transforms.h System.h numerical.h CAP.h molden_transform.h read_qchem_fchk.h
OBJ = Shell.o BasisSet.o Atom.o overlap.o libCAP.o utils.o transforms.o System.o numerical.o molden_transform.o read_qchem_fchk.o


%.o: %.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

main: $(OBJ) 
	$(CXX) $(CXXFLAGS) -o main $(OBJ) $(LIBS)

.PHONY: clean

clean:
	rm *.o
	

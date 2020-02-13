CXX=g++
CXXFLAGS= -Wall -g -std=c++17
LIBS=-larmadillo -lm 
DEPS = Shell.h overlap.h Atom.h BasisSet.h utils.h transforms.h
OBJ = Shell.o BasisSet.o Atom.o overlap.o libCAP.o utils.o transforms.o


%.o: %.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

main: $(OBJ) 
	$(CXX) $(CXXFLAGS) -o main $(OBJ) $(LIBS)

.PHONY: clean

clean:
	rm *.o
	

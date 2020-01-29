CXX=g++
CXXFLAGS= -Wall -g -std=c++11
DEPS = BasisFunction.h overlap.h Atom.h BasisSet.h
OBJ = BasisSet.o Atom.o overlap.o BasisFunction.o libCAP.o 
LIBS=-lm


%.o: %.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

main: $(OBJ) 
	$(CXX) $(CXXFLAGS) -o main $(OBJ) 

.PHONY: clean

clean:
	rm *.o
	

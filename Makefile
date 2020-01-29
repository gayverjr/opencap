CXX=g++
CXXFLAGS= -Wall -g -std=c++11
DEPS = BasisFunction.h overlap.h
OBJ = overlap.o BasisFunction.o libCAP.o 
LIBS=-lm


%.o: %.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

main: $(OBJ) 
	$(CXX) $(CXXFLAGS) -o main $(OBJ) 

.PHONY: clean

clean:
	rm *.o
	

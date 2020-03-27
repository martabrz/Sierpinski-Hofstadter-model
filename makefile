CXX = g++
CXXFLAGS =  -O3 -L/usr/local/opt/openblas/lib -lopenblas -llapack -march=native $\
							-std=c++11 -larmadillo
HEADERS =  CreateSierpinskiCarpet.h CreateHamiltonian.h LocalChern.h
SOURCES =  CreateSierpinskiCarpet.cpp CreateHamiltonian.cpp LocalChern.cpp main.cpp

OUTPUT  = main

all: $(OUTPUT)

$(OUTPUT): $(SOURCES) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $(OUTPUT) $(SOURCES)

.PHONY: clean
clean:
	rm -f main
	rm -f xy_positions.dat

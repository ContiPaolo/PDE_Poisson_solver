CXXFLAGS = -g -lumfpack
CXX = c++
PROGS = main
all: $(PROGS)

main: main.o Validation.o
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LIBS)

Validation.o: Validation.hpp
main.o: Validation.hpp Stencils.hpp PlotFunction.hpp

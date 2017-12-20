.PHONY: clean
MATLABDIR ?= //usr/local/MATLAB/R2014a
# It should be the location of where matlab is installed
CXX ?= g++
CFLAGS = -std=c++11 -w -fopenmp -lm -Wall -Wconversion -O3 -fPIC -I$(MATLABDIR)/extern/include -I..
LDFLAGS = -fopenmp
LIBS = blas/blas.a
CC ?= gcc
MEX = $(MATLABDIR)/bin/mex
MEX_OPTION = CC="$(CXX)" CXX="$(CXX)" CFLAGS="$(CFLAGS)" CXXFLAGS="$(CFLAGS)" LDFLAGS="$(LDFLAGS)"
# comment the following line if you use MATLAB on a 32-bit computer
MEX_OPTION += -largeArrayDims
MEX_EXT = $(shell $(MATLABDIR)/bin/mexext)
MAXCUTdir ?= MixingSDPSolve
MAXCUTlocal ?= $(MAXCUTdir)/MixMaxCutSparseAAT.cpp
all: clean w_solver_recom.$(MEX_EXT) MixMaxCutSparseAAT.$(MEX_EXT)

w_solver_recom.$(MEX_EXT): w_solver_recom.cpp tron.o blas/blas.a
	$(MEX) $(MEX_OPTION) w_solver_recom.cpp tron.o blas/blas.a
MixMaxCutSparseAAT.$(MEX_EXT): $(MAXCUTlocal)
	$(MEX) $(MEX_OPTION) $(MAXCUTlocal) -I$(MAXCUTdir) 
tron.o: tron.cpp tron.h
	$(CXX) $(CFLAGS) -c -o tron.o tron.cpp

blas/blas.a: blas/*.c blas/*.h
	make -C blas OPTFLAGS='$(CFLAGS)' CC='$(CC)';

clean:
	make -C blas clean
	rm -f *~ *.o *.$(MEX_EXT)

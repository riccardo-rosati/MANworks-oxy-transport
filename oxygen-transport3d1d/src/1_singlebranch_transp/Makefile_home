# ====================================================================
#   "Mixed Finite Element Methods for Coupled 3D/1D Transport Problems"
#      Course on Advanced Programming for Scientific Computing
#                     Politecnico di Milano
#                         A.Y. 2015-2016
#
#                    Copyright (C) 2016 Stefano Brambilla
# ====================================================================
#   FILE        : Makefile
#   DESCRIPTION : makefile for test simulations
#   AUTHOR      : Stefano Brambilla <s.brambilla93@gmail.com> 
#   DATE        : September 2016
# ====================================================================

#Da utilizzare su macchina virtuale, altrimenti usare Makefile_MOX

include ../../config.mk

CPPFLAGS=-I../../include -I$(GETFEM_PREFIX)/include -I$(PROBLEM_FLUID_PREFIX)/include -I$(mkBoostInc) -DGMM_USES_LAPACK -DGMM_USES_BLAS
CXXFLAGS=-std=c++14 -D=M3D1D_VERBOSE_

ifeq ($(DEBUG),yes)
  OPTFLAGS=-g -Wall
else
  OPTFLAGS=-O3 -march=native
  CPPFLAGS+=-DNDEBUG
endif

LDFLAGS+=-L../../lib -L$(GETFEM_PREFIX)/lib -L$(PROBLEM_FLUID_PREFIX)/lib -Wl,-rpath=../../lib -Wl,-rpath=$(PROBLEM_FLUID_PREFIX)/lib -L$(mkBoostLib)
LDFLAGS+=-L../../utilities/muparser/lib -Wl,-rpath=../../utilities/muparser/lib 
#LDFLAGS += -L/opt/lib/samg/
LIBRARIES=-lproblem3d1d -lgetfem -loxygen_transport3d1d -lutil -lboost_iostreams -lboost_system -lboost_filesystem -L/home/riccardo/Documenti/software/blas/blas/lib -lblas -llapack -L$(mkQhullLib) -lqhull
#LIBRARIES += -lamg -liomp5 

SRCS=$(wildcard *.cpp)
OBJS=$(SRCS:.cpp=.o)
EXEC=M3D1D


OUTDIR=vtk

.PHONY: all clean distclean

all: $(EXEC)
	@echo
	@echo Compilation completed!

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPTFLAGS) -o $@ -c $<

$(EXEC): $(OBJS)
	$(CXX) $(OBJS) -o $@  $(LDFLAGS) $(LIBRARIES)

clean:
	$(RM) $(OBJS) $(EXEC) *~ *.log

distclean: clean
	$(RM) *.txt $(OUTDIR)/*

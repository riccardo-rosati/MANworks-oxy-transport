# Mixed Finite Element Methods for Coupled 3D/1D Fluid Problems
#### *Politecnico di Milano* (ITALY)

**Authors** :  Luca Possenti, Simone Di Gregorio, Giorgio Raimondi, Fannie Gerosa

**Mailto** : <luca.possenti@polimi.it>

**Date**   : January 2018

-------------------------------------------------------
## How to install and run the program
-------------------------------------------------------
## THE PACKAGE

- `include/` : General include files

- `lib/`     : Main library (to be generated)

- `src/`     : Example sources
  - 2_singlebranch/ : solve the coupling with single-vessel network
  - 3_Ybifurcation/ : solve the problem with Y-shaped bifurcation
  - 7_curved_singleBranch: solve the problem with a curved single-vessel network
  - 8_curved_bifurcation: solve the problem with Y-shaped bifurcation
  - 9_voronoi_network: solve the problem with a complex network

- `config.mk`: Specify the variable GETFEM_PREFIX for GetFEM++ library

- `Makefile` : Instruction to install the project (see INSTALLATION)

- `utilities`: other files for problem set up

## INSTALLATION
### Prerequisites

You need the open source finite element library "GetFEM++"

See <http://download.gna.org/getfem/html/homepage>

Version >= 5.0 is preferible

You have to modify the path to the GetFEM library in `config.mk`:
``` 
GETFEM_PREFIX=/home/.../path/to/.../getfem
``` 

Alternatively, at MOX cluster use the `module.sh` file:
``` 
$ source configure.sh
``` 
Also SAMG LECENCE is required.

Gnuplot: 
Gnuplot is NOT required, but it can be used to visualize residuals. 
To use it, uncomment lines within the code and see the GNUPLOT_Istruzioni_installazione to install.
https://sourceforge.net/projects/gnuplot/files/gnuplot/

BEWARE: 
Recall to add the library path to LD_LIBRARY_PATH. Example:
```
$ export LD_LIBRARY_PATH=/home/...path/to.../getfem/lib

```
======================

### Installation
Build the whole project with:
``` 
$ make
``` 
It first build the (static) library "libproblem3d1d" by calling
the Makefile in `include/`:
``` 
$ make -C include/
``` 
Then, it calls the inner makefiles provided for all examples.

It is also possible to build a single example, e.g. "2_singlebranch", with:
``` 
$ make library

$ make -C src/2_singlebranch
``` 

BEWARE: 
If you want non-optimized program type:
``` 
$ make DEBUG=yes 
``` 
By defaul DEBUG=no.

The following macro are defined and exported
``` 
CPPFLAGS=-I../../include -I$(GETFEM_PREFIX)/include

CXXFLAGS=-std=c++11 -D=M3D1D_VERBOSE_

OPTFLAGS=-O3 -DNDEBUG -march=native

LDFLAGS=-L$(GETFEM_PREFIX)/lib

LIBRARIES=-lgetfem
``` 
Recall that any macro may be overrlued by specifying it when calling 
make. Example: 
``` 
$ make CXXFLAGS+=-DSOMETHING OPTFLAGS=-g
``` 

======================

### Documentation (only implemented for Notaro's code)
The documentation is produced by doxygen. The file Doxyfile contains 
the common doxygen configuration for all examples.
Build the code documentation with:
``` 
$ make pdf
``` 
It first fills doc/ with code documentation ($ make doc) and then compile
the .tex files to produce a portable file ($ pdflatex doc/latex/refman.tex).
You can visualize the documentation with
``` 
$ evince doc/latex/refman.pdf
``` 

## MAKE OPTIONS
All examples are provided with a Makefile which accepts the following
options:
-  all       : makes the example
-  clean     : as it says
-  distclean : clean and also deletes temporary file and local doc directory
Being "all" the first target of the makefile, to compile the examples is
sufficient to type make. 
In addition the external Makefile (./Makefile) has the following options:
-  doc       : produces the documentation (html, tex)
-  pdf       : produces a guide in portable format
- library    : build the library from files in include/

## RUN EXAMPLES
To run a specific example, go to the related subdirectory
``` 
$ cd src/2_singlebranch
``` 
Build the program
``` 
$ make
``` 
Execute the program with specific input
``` 
$ ./M3D1D input.param
``` 
Each program contains the file input.param that specifies 

- Some flags to identify the particular example
  -  TEST_PARAM = 1  # import parameters in dimensionless form
  -  VTK_EXPORT = 1  # export results in vtk format
  -  ...

- The mesh
  - For the 3D mesh you can either provide instruction to build a simple
  regular mesh (TEST_GEOMETRY = 1) or the absolute path to import a mesh
  pre-built with Gmsh (.gmsh)
  - For the 1D mesh specify the path to the file of points (.pts). All
  examples come with a possible pts file

- GetFEM++ descriptors (FEM, ...)

- Problem parameters (dimensional or dimensionless)

- Boundary conditions. You can choose an arbitrary combination of
  Dirichlet-type conditions on pt and/or Robin-type conditions
  on the flux, namely:

  % Faces:   x=0  x=L  y=0  y=L  z=0  z=L

  % BC labels (DIR / MIX)

  BClabel = 'DIR  DIR  DIR  DIR  DIR  DIR'

  % BC values

  BCvalue = '0.0  0.0  0.0  0.0  0.0  0.0'
  

##  DEV ENVIRONMENT
OS         : Ubuntu 14.04 LTS 64-bit (Tested also on Centos)

Compiler   : g++-5.2.1

GetFEM lib : 5.0

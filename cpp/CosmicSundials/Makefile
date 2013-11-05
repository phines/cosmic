#file: makefile for Cosmic
#author: 2012 Ben O'Hara
#contact: buohara@uvm.edu, buohara@gmail.com
#description: makefile for the main Cosmic program

#pull in the configuration file
include config.mk

#default target is the serial debug version (for now)
default: CosmicDS

FILES = kin_klu psc utils ic_solver_NR
OBJS = $(addsuffix .o, $(addprefix lib/, $(FILES)))
INC = $(addsuffix .h, $(addprefix include/, $(FILES)))

I = -Iinclude
LIBS = -lklu -lamd -lcolamd -lbtf -lsundials_kinsol -lcxsparse
LIBS_S = $(LIBS) -lsundials_nvecserial
LIBS_P = $(LIBS) -lsundials_nvecparallel

#short-hand for the various sets of compiler options
LFDS = $(LFLAGS_D) $(I) $(OBJS) $(LIBS_S)
LFS = $(LFLAGS) $(I) $(OBJS) $(LIBS_S)
LFDP = $(LFLAGS_D) $(MPI) $(I) $(OBJS) $(LIBS_P)
LFP = $(LFLAGS) $(MPI) $(I) $(OBJS) $(LIBS_P)

#debug serial. links the serial nvector module
CosmicDS: library $(INC)
	$(CC) $(LFDS) src/main.c -o CosmicDS
	
#release serial. links the nvector serial module and turns on optimizations
CosmicS: library $(INC)
	$(CC) $(LFS) src/main.c -o CosmicS

#debug parallel. defines MPI and links the nvector parallel module
CosmicDP: library $(INC)
	$(CC) $(LFDP) src/main.c -o CosmicDP

#release parallel. defines MPI, links the nvector parallel module, and turns on optimizations
CosmicP: library $(INC)
	$(CC) $(LFP) src/main.c -o CosmicP

#build library files created by the cosmic module	
library:
	(cd lib; $(MAKE))
	
#sandbox testing
sandbox: library $(INC) src/sandbox.c
	gcc -Wall $(OBJS) src/sandbox.c $(LIBS_S) -o Sandbox
	
#remove cosmic libraries and executables
clean:
	(cd lib; $(MAKE) clean)
	rm -f Cosmic* Sandbox
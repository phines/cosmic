#file: Makefile for CosmicSundials libraries
#author: 2012 Ben O'Hara
#contact: buohara@gmail.com, buohara@uvm.edu
#description: makefile for building various cosmic-related libraries.

include ../config.mk

I = -I../include
VPATH = ../src:../include
OBJS = kin_klu.o psc.o utils.o ic_solver_NR.o

default: all

all: $(OBJS)

ic_solver.o: ic_solver.c ic_solver.h psc.o kin_klu.o
	$(CC) $(CFLAGS_D) $(I) ../src/ic_solver.c
	
ic_solver_NR.o: ic_solver_NR.c ic_solver_NR.h psc.o
	$(CC) $(CFLAGS_D) $(I) ../src/ic_solver_NR.c
	
psc.o: psc.c psc.h
	$(CC) $(CFLAGS_D) $(I) ../src/psc.c
	
kin_klu.o: kin_klu.c
	$(CC) $(CFLAGS_D) $(I) ../src/kin_klu.c
	
%.o: %.c %.h
	$(CC) $(CFLAGS_D) $(I) ../src/$< -o $@
	
clean:
	rm -f *.o
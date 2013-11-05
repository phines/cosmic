/**
 file: main.cpp
 author: Ben O'Hara
 email: buohara@gmail.com
 description: main entry point for the cosmic cascading failure simulator.
 */

#include <stdio.h>
#include <stdlib.h>
#include "../include/psc.h"
#include "../include/ic_solver_NR.h"

int main(int argc, char **argv){

	char* file_name;
	PS ps;
	double *xy0;

	printf("Launching Cosmic...\n");

	//initialize MPI if we are using parallel code.
#ifdef HAVE_MPI
	MPI_Init(&argc, &argv);
#endif
	
	/*
	get the ps file from the terminal. assume the first argument is
	always the ps file.
	 */
	if(argc>1)
		file_name=argv[1];
	else{
		printf("Please specify PS data file.");
		return 1;
	}
	
	//get PS data from the file
	ps=GetPS(file_name);
	
	//check whether we were able to open the file
	if(ps==NULL){
		printf("Unable to build PS data.\n");
		return 1;
	}
	xy0=malloc(sizeof(int));

	//call the initial condition solver.
	SolveIC(ps, xy0);
	
	//free initial condition memory
	free(xy0);

	//terminate MPI
#ifdef HAVE_MPI
	MPI_Finalize();
#endif
	
	//free power system data and return
	FreePS(&ps);
	return 0;
}

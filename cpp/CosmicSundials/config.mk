#file: config.mk
#author: 2012 Ben O'Hara
#contact: buohara@gmail.com, buohara@uvm.edu
#description: make configurations for Cosmic

#compiler settings
CC = gcc
CFLAGS_D = -Wall -c -D_VERBOSE
LFLAGS_D = -Wall -D_VERBOSE
CFLAGS = -Wall -O2 -c
LFLAGS = -Wall -O2
MPI = -DHAVE_MPI

#install directories
INSTALL_LIB = /usr/local/lib
INSTALL_INCLUDE = /usr/local/include
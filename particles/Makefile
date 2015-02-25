#
# Hopper - NERSC 
#
# Portland Group Compilers PGI are loaded by default; for other compilers please check the module list
#

# This makefile is messy.  I apologize for having no idea what I'm doing.  This
# all seems needlessly complicated in comparison to Java build tools.

# Identify what machine we're on and set up some environment-specific flags:
ifneq ($(shell uname -a | grep Ubuntu),)
  RAW_CC = gcc
  OSFOUND=1
endif
ifneq ($(shell uname -a | grep Darwin),)
  RAW_CC = g++-4.9 -Wall -std=c++0x
	OPENMP = -fopenmp
	RAW_MPCC = OMPI_MPICXX=$(RAW_CC) mpic++
  FLAGS = -DCUSTOM_PTHREAD_BARRIER=1
	OWN_PTHREAD_BARRIER_LIB = pthread_barrier.o
  OSFOUND=1
endif
ifneq ($(shell uname -a | grep hopper),)
  HOPPER = 1
  RAW_CC = cc -Wall
	OPENMP = -mp
	RAW_MPCC = $(RAW_CC)
  OSFOUND=1
endif

CC = $(RAW_CC) $(FLAGS)
MPCC = $(RAW_MPCC) $(FLAGS)
ifdef DEBUG
	CFLAGS = -O0 -g -pg
else
	CFLAGS = -O3
endif
LIBS = -lstdc++

OWN_LIBS = common.o Grid.o Stats.o $(OWN_PTHREAD_BARRIER_LIB)

TARGETS = serial pthreads autograder tests
#openmp mpi

all:	$(TARGETS)

tests: tests.o Grid.o
	$(CC) -o $@ tests.o Grid.o $(LIBS)
serial: serial.o $(OWN_LIBS)
	$(CC) -o $@ $^ $(LIBS)
autograder: autograder.o $(OWN_LIBS)
	$(CC) -o $@ $^ $(LIBS)
pthreads: pthreads.o $(OWN_LIBS)
	$(CC) -o $@ $^ $(LIBS) -lpthread
# openmp: openmp.o $(OWN_LIBS)
# 	$(CC) -o $@ $^ $(LIBS) $(OPENMP)
# mpi: mpi.o common.o
# 	$(MPCC) -o $@ $^ $(LIBS) $(MPILIBS)

pthread_barrier.o: pthread_barrier.cpp pthread_barrier.h
	$(CC) -c $(CFLAGS) $<
autograder.o: autograder.cpp $(OWN_LIBS)
	$(CC) -c $(CFLAGS) $<
# openmp.o: openmp.cpp $(OWN_LIBS)
# 	$(CC) -c $(OPENMP) $(CFLAGS) $<
serial.o: serial.cpp $(OWN_LIBS)
	$(CC) -c $(CFLAGS) $<
pthreads.o: pthreads.cpp $(OWN_LIBS)
	$(CC) -c $(CFLAGS) $<
# mpi.o: mpi.cpp $(OWN_LIBS)
# 	$(MPCC) -c $(CFLAGS) $<
common.o: common.cpp common.h Stats.h
	$(CC) -c $(CFLAGS) $<
Stats.o: Stats.cpp Stats.h
	$(CC) -c $(CFLAGS) $<
tests.o: tests.cpp Grid_test.cpp Grid.h catch.hpp
	$(CC) -c $(CFLAGS) $<
%.o: %.cpp %.h
	$(CC) -c $(CFLAGS) $<

clean:
	rm -f *.o $(TARGETS) *.stdout *.txt

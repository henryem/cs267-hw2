#
# Hopper - NERSC 
#
# Portland Group Compilers PGI are loaded by default; for other compilers please check the module list
#

ifneq ($(shell uname -a | grep Ubuntu),)
  RAW_CC = gcc
  OSFOUND=1
endif
ifneq ($(shell uname -a | grep Darwin),)
  RAW_CC = g++-4.9
	OPENMP = -fopenmp
	RAW_MPCC = OMPI_MPICXX=$(RAW_CC) mpic++
  FLAGS = -DCUSTOM_PTHREAD_BARRIER=1
	OWN_PTHREAD_BARRIER_LIB = pthread_barrier.o
  OSFOUND=1
endif
ifneq ($(shell uname -a | grep hopper),)
  HOPPER = 1
  RAW_CC = cc
	OPENMP = -mp
	RAW_MPCC = $(RAW_CC)
  OSFOUND=1
endif

CC = $(RAW_CC) $(FLAGS)
MPCC = $(RAW_MPCC) $(FLAGS)
CFLAGS = -O3
LIBS =

OWN_LIBS = common.o $(OWN_PTHREAD_BARRIER_LIB)

TARGETS = serial pthreads openmp mpi autograder

all:	$(TARGETS)

serial: serial.o $(OWN_LIBS)
	$(CC) -o $@ $(LIBS) serial.o $(OWN_LIBS)
autograder: autograder.o $(OWN_LIBS)
	$(CC) -o $@ $(LIBS) autograder.o $(OWN_LIBS)
pthreads: pthreads.o $(OWN_LIBS)
	$(CC) -o $@ $(LIBS) -lpthread pthreads.o $(OWN_LIBS)
openmp: openmp.o common.o
	$(CC) -o $@ $(LIBS) $(OPENMP) openmp.o $(OWN_LIBS)
mpi: mpi.o common.o
	$(MPCC) -o $@ $(LIBS) $(MPILIBS) mpi.o common.o

pthread_barrier.o: pthread_barrier.cpp pthread_barrier.h
	$(CC) -c $(CFLAGS) pthread_barrier.cpp
autograder.o: autograder.cpp $(OWN_LIBS)
	$(CC) -c $(CFLAGS) autograder.cpp
openmp.o: openmp.cpp $(OWN_LIBS)
	$(CC) -c $(OPENMP) $(CFLAGS) openmp.cpp
serial.o: serial.cpp $(OWN_LIBS)
	$(CC) -c $(CFLAGS) serial.cpp
pthreads.o: pthreads.cpp $(OWN_LIBS)
	$(CC) -c $(CFLAGS) pthreads.cpp
mpi.o: mpi.cpp $(OWN_LIBS)
	$(MPCC) -c $(CFLAGS) mpi.cpp
common.o: common.cpp
	$(CC) -c $(CFLAGS) common.cpp

clean:
	rm -f *.o $(TARGETS) *.stdout *.txt

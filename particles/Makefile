# This makefile is messy.  I apologize for having no idea what I'm doing.  This
# all seems needlessly complicated in comparison to Java build tools.

# Identify what machine we're on and set up some environment-specific flags:
ifneq ($(shell uname -a | grep Ubuntu),)
  # Not actually supported.
  RAW_CC = gcc
  OSFOUND=1
endif
ifneq ($(shell uname -a | grep Darwin),)
  RAW_CC = g++-4.9 -Wall -std=c++0x
  OPENMP = -fopenmp
  RAW_MPCC = OMPI_MPICXX="$(RAW_CC)" mpic++
  FLAGS = -DCUSTOM_PTHREAD_BARRIER=1
  OWN_PTHREAD_BARRIER_LIB = pthread_barrier.o
  OSFOUND=1
  RAW_NVCCFLAGS = -arch=sm_30 # My Mac supports Compute Level 3.0.  Not sure how to find this automatically.
  SUPPORTS_GPU=1
endif
ifneq ($(shell uname -a | grep hopper),)
  # Does not support GPU build.  Also, the build assumes that GNU g++ is the compiler.
  # Run ``module swap PrgEnv-pgi PrgEnv-gnu'' to switch.
  HOPPER = 1
  RAW_CC = cc -std=c++0x -Wall
  OPENMP = -fopenmp
  RAW_MPCC = $(RAW_CC)
  OSFOUND=1
endif
ifneq ($(shell uname -a | grep edison),)
  # Does not support GPU build.  Also, the build assumes that GNU g++ is the compiler.
  # Run ``module swap PrgEnv-pgi PrgEnv-gnu'' to switch.
  HOPPER = 1
  RAW_CC = cc -std=c++0x -Wall
  OPENMP = -fopenmp
  RAW_MPCC = $(RAW_CC)
  OSFOUND=1
endif
ifneq ($(shell uname -a | grep stampede),)
  # Note: Need to ``module load cuda'' first.
  STAMPEDE = 1
  RAW_CC = cc -Wall
  OPENMP = -mp
  RAW_MPCC = $(RAW_CC)
  OSFOUND=1
  RAW_NVCCFLAGS = -arch=sm_35
  SUPPORTS_GPU=1
endif

CC = $(RAW_CC) $(FLAGS)
MPCC = $(RAW_MPCC) $(FLAGS)
ifdef DEBUG
  CFLAGS = -O0 -g -pg
  RAW_NVCCFLAGS += -DTHRUSTDEBUG
else
  CFLAGS = -O3
endif
LIBS = -lstdc++

# NVCC stuff:
NVCC = nvcc
NVCCFLAGS = $(CFLAGS) $(RAW_NVCCFLAGS)
NVCC_LIBS = -lm

OWN_LIBS = common.o Grid.o Stats.o $(OWN_PTHREAD_BARRIER_LIB)
OWN_OMP_LIBS = $(OWN_LIBS) OmpThreadsafeGrid.o

TARGETS = serial serial-naive autograder tests openmp mpi
ifdef SUPPORTS_GPU
  TARGETS += gpu gpu-naive
endif

all:	$(TARGETS)

tests: tests.o Grid.o
	$(CC) $(CFLAGS) -o $@ tests.o Grid.o $(LIBS)
serial: serial.o $(OWN_LIBS)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)
serial-naive: serial-naive.o $(OWN_LIBS)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)
autograder: autograder.o $(OWN_LIBS)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)
pthreads: pthreads.o $(OWN_LIBS)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS) -lpthread
openmp: openmp.o $(OWN_OMP_LIBS)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS) $(OPENMP)
mpi: mpi.o $(OWN_LIBS)
	$(MPCC) $(CFLAGS) -o $@ $^ $(LIBS) $(MPILIBS)
gpu: gpu.o common_gpu.o Stats.o
	$(NVCC) -o $@ $(NVCCLIBS) gpu.o common_gpu.o Stats.o
gpu-naive: gpu_naive.o common_gpu.o Stats.o
	$(NVCC) -o $@ $(NVCCLIBS) gpu_naive.o common_gpu.o Stats.o

# Load CUDA using the following command
# module load cuda

pthread_barrier.o: pthread_barrier.cpp pthread_barrier.h
	$(CC) -c $(CFLAGS) $<
autograder.o: autograder.cpp $(OWN_LIBS)
	$(CC) -c $(CFLAGS) $<
openmp.o: openmp.cpp $(OWN_OMP_LIBS)
	$(CC) -c $(OPENMP) $(CFLAGS) $<
serial.o: serial.cpp $(OWN_LIBS)
	$(CC) -c $(CFLAGS) $<
serial-naive.o: serial-naive.cpp $(OWN_LIBS)
	$(CC) -c $(CFLAGS) $<
pthreads.o: pthreads.cpp $(OWN_LIBS)
	$(CC) -c $(CFLAGS) $<
mpi.o: mpi.cpp $(OWN_LIBS)
	$(MPCC) -c $(CFLAGS) $<
common.o: common.cpp common.h Stats.h
	$(CC) -c $(CFLAGS) $<
OmpThreadsafeGrid.o: OmpThreadsafeGrid.cpp OmpThreadsafeGrid.h
	$(CC) -c $(OPENMP) $(CFLAGS) $<
tests.o: tests.cpp Grid_test.cpp Grid.h catch.hpp
	$(CC) -c $(CFLAGS) $<
gpu.o: gpu.cu common.h
	$(NVCC) -c $(NVCCFLAGS) gpu.cu
gpu_naive.o: gpu_naive.cu common.h
	$(NVCC) -c $(NVCCFLAGS) gpu_naive.cu
common_gpu.o: common.cpp common.h
	$(NVCC) --x=cu -o $@ -c $(NVCCFLAGS) common.cpp
%.o: %.cpp %.h
	$(CC) -c $(CFLAGS) $<

clean:
	rm -f *.o $(TARGETS) *.stdout *.txt

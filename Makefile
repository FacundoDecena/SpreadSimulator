# Build Executable
SHELL=/bin/bash
# constants
np = 16
npl = 6
mph = 2
CC = gcc
MPCC = mpicc
MPR = mpirun

clean:
	$(RM) ./build/*
	touch ./build/.nofiles

build:
	${CC} 	-o 	./build/secuencial	./secuencial/main.c -fopenmp
	${CC} 	-o 	./build/omp 		./omp/main.c 		-fopenmp
	${MPCC}	-o 	./build/mpi 		./mpi/main.c  		-fopenmp
	${MPCC}	-o 	./build/hybrid 		./hybrid/main.c  	-fopenmp

.PHONY: clean build

secuencial:
	./build/secuencial

omp:
	./build/omp

mpi:
	${MPR} -np ${npl} ./build/mpi
	
hybrid:
	${MPR} -np ${mph} ./build/hybrid

.PHONY: secuencial omp mpi hybrid
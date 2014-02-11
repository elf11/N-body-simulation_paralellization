CFLAGS=-g -Wall

all: serial omp pthread mpi

serial: nbody_serial.cpp
	g++ $(CFLAGS) -o serial nbody_serial.cpp -fopenmp

omp: nbody_omp.cpp
	g++ $(CFLAGS) -o omp nbody_omp.cpp -fopenmp

pthread: nbody_threads.cpp
	g++ $(CFLAGS) -o pthread nbody_threads.cpp -fopenmp
	
mpi: nbody_mpi.cpp
	mpic++ $(CFLAGS) -o mpi nbody_mpi.cpp -fopenmp

clean:
	rm -rf serial omp pthread mpi *.out

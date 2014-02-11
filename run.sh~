#!/bin/bash

rm time_values.txt
touch time_values.txt

./serial
./omp
./pthread
mpirun -np 4 ./mpi

cat time_values.txt

#! /bin/bash

# Compile the serial version
gcc -o recover_pins_serial serial.c -lcrypto -lm

# Compile the OpenMP version
gcc -o recover_pins_openmp -fopenmp openmp.c -lcrypto

# Compile the cilkplus version
gcc -o recover_pins_cilk -fcilkplus cilk.c -lcrypto

# Compile the tbb version
g++ -o recover_pins_tbb tbb.c -ltbb

# Compile the MPI version
mpicc -o recover_pins_mpi mpi.c

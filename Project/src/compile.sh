#! /bin/bash

srun g++ -o stencil_serial1 filter1_serial.cpp -lopencv_core -lopencv_imgcodecs
srun g++ -o stencil_openmp1 -fopenmp filter1_omp.cpp -lopencv_core -lopencv_imgcodecs



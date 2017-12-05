srun g++ -o stencil_serial filter2_serial.cpp -lopencv_core -lopencv_imgcodecs  
srun g++ -o stencil_openmp -fopenmp filter2_omp.cpp -lopencv_core -lopencv_imgcodecs



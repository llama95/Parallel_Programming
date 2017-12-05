Our image is from NASA at:
http://cs.coloradocollege.edu/~dellsworth/image.jpg

You can use wget or curl to download these images.


The following modules should be loaded with “module load ____” before working with this code:
slurm
gcc-7.2.0
opencv-3.3.0
openmpi-3.0.0

Use “compile1.sh” to compile the first filter and “compilefilter2.sh” to compile the second filter. 
Enter “./stencil_serial1 image.jpg” to see the modified image with the first filter and “./stencil_serial image.jpg” to see the image modified with the second filter. 
Enter “./stencil_openmp1 image.jpg” to see the image modified in parallel with the first filter and “./stencil_openmp image.jpg” to see the image modified in parallel with the second filter.

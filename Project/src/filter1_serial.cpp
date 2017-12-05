#include<iostream>
#include<cmath>
#include <opencv2/opencv.hpp>
#include<opencv2/imgproc/imgproc.hpp>
#include<opencv2/highgui/highgui.hpp>

using namespace std;
using namespace cv;

struct pixel {
    double red; // class and constructor
    double green;
    double blue;

    pixel(double r, double g, double b) : red(r), green(g), blue(b) {};
};


void filter1x_kernel(const int rows, const int cols, double * const kernel) {
    if(rows != 3 || cols !=3) {
        std::cerr << "Bad Prewitt kernel matrix\n"; //builds a matrix to build the image within the kernal
        return;
    }
    for(int i=0;i<3;i++) {
        kernel[0 + (i*rows)] = -1.0;
        kernel[1 + (i*rows)] = 0.0;
        kernel[2 + (i*rows)] = 1.0;
    }
    kernel[1] = -2;
    kernel[7] = 2;



}

void filter1y_kernel(const int rows, const int cols, double * const kernel) {
        if(rows != 3 || cols !=3) {
                std::cerr << "Bad Prewitt kernel matrix\n";
                return;
        }
        for(int i=0;i<3;i++) {
                kernel[i + (0*rows)] = 1.0;
                kernel[i + (1*rows)] = 0.0;
                kernel[i + (2*rows)] = -1.0;



        }
        kernel[3] = 2;
        kernel[5] = -2;

  }


  void apply_stencil(const int radius, const double stddev, const int rows, const int cols, pixel * const in, pixel * const out) { // how blurry how big image coming in image exiting program as array of pixels
      const int dim = radius*2+1;
      double kernel[dim*dim]; //kernel arrays // take in image and takes in kernel
      filter1x_kernel(dim, dim, kernel);
      filter1y_kernel(dim, dim, kernel);
      // For each pixel in the image...
      for(int i = 0; i < rows; ++i) {
          for(int j = 0; j < cols; ++j) {
              const int out_offset = i + (j*rows);
              // ...apply the template centered on the pixel...
              for(int x = i - radius, kx = 0; x <= i + radius; ++x, ++kx) { // go and look at the pixels around the image given by the kernal
                  for(int y = j - radius, ky = 0; y <= j + radius; ++y, ++ky) {
                      // ...and skip parts of the template outside of the image
                      if(x >= 0 && x < rows && y >= 0 && y < cols) {
                          // Acculate intensities in the output pixel
                          const int in_offset = x + (y*rows);
                          const int k_offset = kx + (ky*dim);
                          out[out_offset].red   += kernel[k_offset] * in[in_offset].red;
                          out[out_offset].green += kernel[k_offset] * in[in_offset].green; // accumulate the color values for each within the pixel
                          out[out_offset].blue  += kernel[k_offset] * in[in_offset].blue;
                      }
                  }
              }
          }
      }
  }

  int main( int argc, char* argv[] ) {
      if(argc != 2) {
          std::cerr << "Usage: " << argv[0] << " imageName\n";
          return 1;
      }

      // Read image
      Mat image;
      image = imread(argv[1],CV_LOAD_IMAGE_COLOR);
      if(!image.data ) {
          std::cout <<  "Error opening " << argv[1] << std::endl;
          return -1;
      }

      // Get image into C array of doubles for processing
      const int rows = image.rows;
      const int cols = image.cols;
      pixel * imagePixels = (pixel *) malloc(rows * cols * sizeof(pixel));
      for(int i = 0; i < rows; ++i) {
          for(int j = 0; j < cols; ++j) {
              Vec3b p = image.at<Vec3b>(i, j);
              imagePixels[i + (j*rows)] = pixel(p[0]/255.0,p[1]/255.0,p[2]/255.0);
          }
      }

      // Create output array
      pixel * outPixels = (pixel *) malloc(rows * cols * sizeof(pixel));
      for(int i = 0; i < rows * cols; ++i) {
          outPixels[i].red = 0.0;
          outPixels[i].green = 0.0;
          outPixels[i].blue = 0.0;
      }

      // Do the stencil
      struct timespec start_time;
      struct timespec end_time;
      clock_gettime(CLOCK_MONOTONIC,&start_time);
      apply_stencil(1, 32.0, rows, cols, imagePixels, outPixels);
      clock_gettime(CLOCK_MONOTONIC,&end_time);
      long msec = (end_time.tv_sec - start_time.tv_sec)*1000 + (end_time.tv_nsec - start_time.tv_nsec)/1000000;
      printf("Stencil application took %dms\n",msec);

      // Create an output image (same size as input)
      Mat dest(rows, cols, CV_8UC3);
      // Copy C array back into image for output
      for(int i = 0; i < rows; ++i) {
          for(int j = 0; j < cols; ++j) {
              const size_t offset = i + (j*rows);
              dest.at<Vec3b>(i, j) = Vec3b(floor(outPixels[offset].red * 255.0),
                                           floor(outPixels[offset].green * 255.0),
                                           floor(outPixels[offset].blue * 255.0));
          }
      }

      imwrite("out.jpg", dest);


      free(imagePixels);
      free(outPixels);
      return 0;
  }


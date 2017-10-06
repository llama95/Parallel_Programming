/** YOOOOO
 * Program to scale and project data into 2D and find a centroid
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <cilk/cilk.h>
#include <cilk/reducer_opadd.h>

/* struct to hold objects attributes */
struct volume { // struc of array
    float* x; // x coords for all of phasballs
    float* y;
    float* z;
    float* mass;
  //int size = 100;
  //int last = 0;
};

/*
struct volume {
    size_t size;
    size_t last;
    struct phaseball** objects;
};

// Add phaseball to a volume
void volume_append(struct volume* v, struct phaseball* o) {
    if( v->last == v->size ) {
        (v->size) += 100;
        v->objects = realloc(v->objects, sizeof(struct phaseball*)*(v->size)+100);
    }
    (v->objects)[(v->last)] = o;
    (v->last) += 1;
    return;
}
*/

// Place phaseballs uniformly in a box, with mass equal to the manhattan distanced
int place_uniformly(int sx, int ex, int sy, int ey, int sz, int ez, struct volume* v) {
   int counter = 0;
    //need a counter for the # of phaseballs
    // how to calc # of phaseballs??
    // need to allocate memory cuz these are now arrays
    // array will point to x,y,z malloc with sizeof float
  int number_of_phaseballs = (ex-sx+1)*(ey-sy+1)*(ez-sz+1);

  // v->size=100;
  //v->last=0;
  // v.objects = malloc(sizeof(struct spaceball*)*100); // bounds of
  v->x = (float*) malloc(number_of_phaseballs * sizeof(float));
  v->y = (float*) malloc(number_of_phaseballs * sizeof(float));
  v->z = (float*) malloc(number_of_phaseballs * sizeof(float));
  v->mass = (float*) malloc(number_of_phaseballs * sizeof(float));
  for(int i=sx; i<=ex; i++) {
        for(int j=sy; j<=ey; j++) {
            for(int k=sz; k<=ez; k++) {
                // struct volume* n = malloc(sizeof(struct volume));
                v->x[counter] = i; //array pts to each coord counter
                v->y[counter] = j; // array will point to x,y,z malloc with sizeof float
                v->z[counter] = k;
                //v->mass[i] = 1; //
                v->mass[counter] = fabs(v->x[counter])+fabs(v->y[counter])+fabs(v->z[counter]);
                // volume_append(v,n); // add to counter
		 counter++;
            }
        }
  }
  return number_of_phaseballs;
}


// Projects 3D volume to 11x11 2D map and report centroid
void post_process(struct volume* v, float* cx, float* cy,int length
) {

  printf("%d\n",length);
  //double mass_sum=0.0;
    //double wx=0.0;
    //double wy=0.0;
    CILK_C_REDUCER_OPADD(mass_sum, double, 0.0);
    CILK_C_REDUCER_OPADD(wx, double, 0.0);
    CILK_C_REDUCER_OPADD(wy, double, 0.0);


    cilk_for(int i=0; i<length; i++) {
        // struct volume* o = v- >objects[i];
      REDUCER_VIEW(mass_sum) += v->mass[i];
      REDUCER_VIEW(wx) += v->x[i] * v->mass[i];
      REDUCER_VIEW(wy) +=v->y[i] * v->mass[i];
    }
    *cx = (wx.value)/(mass_sum.value); //
    *cy = (wy.value)/(mass_sum.value);

    return;
}

int main(int argc, char** argv) {
    // make a volume to store objects in
    // first call should be to our array struct
    struct volume v;
    // Set the initial configuration
    int lengtharray = place_uniformly(-1000,1000,-100,100,-100,100,&v); // include our struct instea//d of &v this will be length of our array

    // Generate output

    float cx,cy;
    //int ex,sx,ey,sy,ez,sz;
    // int number_of_phaseballs = (ex-sx+1)*(ey-sy+1)*(ez-sz+1);
    //int length = number_f_phaseballs;
    struct timespec start_time;
    struct timespec end_time;
    clock_gettime(CLOCK_MONOTONIC,&start_time);
    post_process(&v, &cx, &cy,lengtharray); // lengh of array and our struct
    clock_gettime(CLOCK_MONOTONIC,&end_time);
    long msec = (end_time.tv_sec - start_time.tv_sec)*1000 + (end_time.tv_nsec - start_time.tv_nsec)/1000000;

    // Print output
    printf("Centroid at: %f,%f\n",cx,cy);
    printf("found in %dms\n",msec);

    return EXIT_SUCCESS;
}

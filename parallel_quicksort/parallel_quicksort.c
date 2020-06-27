/*
   Parallel Quicksort - Copyright by Silvio do Lago Pereira, 2019
    
   Intel(R) Core(TM) i7-5500U CPU @ 2.40GHz
   2 cores (4 logical processors)
   RAM 4,00 GB
   Pelles C, v. 8.00.60, 32 bits
   Windows 10, 64 bits

   IMPORTANT: To correctly compile this program, please activate the option
   "Enable OpenMP extensions", by following (from the menu bar):
   Project -> Project Options -> Compiler -> Options
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include <minmax.h>
#include <signal.h>
#include <time.h>
#include <math.h>
#include <omp.h>


// funcion max
#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

// Simple auxiliary inline functions

inline void swap(float v[], int i, int j) { float x=v[i]; v[i]=v[j]; v[j]=x; }
inline double lg(double x, double b) { return log(x)/log(b); }
inline int best_at(int n) { return max(floor(lg(n,4)-4),1); }
inline int best_cutoff(int n) { return max(n>>best_at(n),1); }


// Hoare partition algorithm (avoids worst case for sorted sequences)

int partition(float v[], int start, int end) {
   float pivot = v[(start + end)/2];
   int cut = end+1;
   start--;
   while( start < cut ) {
      do cut--;   while( v[cut]   > pivot );
      do start++; while( v[start] < pivot );
      if( start < cut ) swap(v, start, cut);
   }
   return cut;
}


// Serial Standard Quicksort (may cause stack overflow)

void ssqs(float v[], int start, int end) {
   if( start >= end ) return;
   int cut = partition(v, start, end);
   ssqs(v, start, cut);
   ssqs(v, cut+1, end);
}

void SSqs(float v[], int n) {
   ssqs(v, 0, n-1);
}


// PSqs - Parallel Standard Quicksort

void psqs(float v[], int start, int end, int cutoff) {
   if( end-start+1 <= cutoff )
      ssqs(v, start, end);
   else {
      int cut = partition(v, start, end);
      #pragma omp task
      psqs(v, start, cut, cutoff);
      #pragma omp task
      psqs(v, cut+1, end, cutoff);
   }
}

void PSqs(float v[], int n, int cutoff) {
   #pragma omp parallel
   #pragma omp single nowait
   psqs(v, 0, n-1, cutoff);
}


// FPSqs - Fine-grained Parallel Standard Quicksort

void FPSqs(float v[], int n) {
   #pragma omp parallel
   #pragma omp single nowait
   psqs(v, 0, n-1, 1);
}


// CPSqs - Coarse-grained Parallel Standard Quicksort

void CPSqs(float v[], int n) {
   #pragma omp parallel
   #pragma omp single nowait
   psqs(v, 0, n-1, n>>1);
}


// BPSqs - Best-grained Parallel Standard Quicksort

void BPSqs(float v[], int n) {
   #pragma omp parallel
   #pragma omp single nowait
   psqs(v, 0, n-1, best_cutoff(n));
}


// SOqs  - Serial Optimized Quicksort

void soqs(float v[], int start, int end) {
   while( start < end ) {
      int cut = partition(v, start, end);
      if( cut-start+1 < end-cut ) {
         soqs(v, start, cut);
         start = cut + 1;
      }
      else {
         soqs(v, cut+1, end);
         end = cut;
      }
   }
}

void SOqs(int v[], int n) {
   soqs(v, 0, n-1);
}


// POqs - Parallel Optimized Quicksort

void poqs(float v[], int start, int end, int cutoff) {
   while( end-start+1 > cutoff ) {
      int cut = partition(v, start, end);
      if( cut-start+1 < end-cut ) {
         #pragma omp task 
         poqs(v, start, cut, cutoff);
         start = cut+1;
      }
      else {
         #pragma omp task
         poqs(v, cut+1, end, cutoff);
         end = cut;
      }
   }
   soqs(v, start, end); 
}

void POqs(float v[], int n, int cutoff) {
   #pragma omp parallel
   #pragma omp single nowait
   poqs(v, 0, n-1, cutoff);
}


// FPOqs - Fine-grained Parallel Optimized Quicksort

void FPOqs(float v[], int n) {
   #pragma omp parallel
   #pragma omp single nowait
   poqs(v, 0, n-1, 1);
}


// CPOqs - Coarse-grained Parallel Optimized Quicksort

void CPOqs(float v[], int n) {
   #pragma omp parallel
   #pragma omp single nowait
   poqs(v, 0, n-1, n>>1);
}


// BPOqs - Best-grained Parallel Optimized Quicksort

void BPOqs(int v[], int n) {
   #pragma omp parallel
   #pragma omp single nowait
   poqs(v, 0, n-1, best_cutoff(n));
}


// CSLqs - C Standard Library Quicksort

int ascending(const void *a, const void *b) {
   return (*(int*)a > *(int*)b)-(*(int*)b > *(int*)a);
}

void CSLqs(int v[], int n) {
   qsort(v, n, sizeof(int), ascending);
}


// Auxiliary functions for testing and debugging

void permutation(int v[], int n) {
   for(int i=0; i<n; i++) v[i] = i+1;
   while( n-- ) swap(v, n, rand()%(n+1));
}

void uniform(int v[], int n) {
   for(int i=0; i<n; i++) 
      v[i] = rand()%n;
}

double etime(void *f, float *v, int n, float *w, int cutoff) {
   memcpy(w, v, n*sizeof(int));
   double start = omp_get_wtime();
   if( cutoff <= 0 ) ((void (*)(int *, int))f)(w, n);
   else ((void (*)(int *, int, int))f)(w, n, cutoff);
   return omp_get_wtime()-start;
}

int sorted(int v[], int n) {
   for(int i=0; i<n-1; i++)
      if( v[i] > v[i+1] ) 
         return 0;
   return 1;
}

void show(int v[], int n) {
   for(int i=0; i<n; i++)
      printf("%d ", v[i]);
   puts("\n");
}

void error(int s) {
   printf("Segmentation fault (%d)\n", s);
   abort();
}


// Experiments

// FIND CUTOFF: find point where the speedup is maximized

void find_cutoff(void (*fill)(int *, int), int rounds) {
   double start_experiment = omp_get_wtime();

   for(int k=10; k<=23; k++) {
      int n = pow(2,k);
      int C = best_at(n);

      printf("# Average-time (and speedup) for sorting 2^%d=%d items, using %d threads and %d rounds (best at c=%d)\n", k, n, omp_get_max_threads(), rounds, C);
      printf("#\n");
      printf("#  c     SSqs     POqs  speedup\n");

      double start = omp_get_wtime();

      int *v = malloc(n*sizeof(int));
      int *w = malloc(n*sizeof(int));

      for(int c=1; c<=k; c++ ) {
         int cutoff = max(n>>c,1);

         double tSSqs = 0;
         double tPOqs = 0;

         for(int r=1; r<=rounds; r++) {   
            fill(v, n);
            tSSqs += etime(SSqs, v, n, w, 0);
            tPOqs += etime(POqs, v, n, w, cutoff);
         }

         tSSqs /= rounds;
         tPOqs /= rounds;

         printf("%4d %8.5f %8.5f %8.2f\n", c, tSSqs, tPOqs, tSSqs/tPOqs);
      }
   
      free(v);
      free(w);

      printf("\n# Total elapsed time: %5.2fs\n\n\n", omp_get_wtime()-start);
   }

   printf("\n# Total elapsed time for FIND-CUTOFF experiment: %5.2fs\n\n\n", omp_get_wtime()-start_experiment);
}


// PERFORMANCE: comparison of all Quicksort's versions

void performance(float *v, int rounds) {
   printf("# Average-time (and speedup) for sorting n=2^k items, using %d threads and %d rounds\n",omp_get_max_threads(),rounds);
   printf("#\n");
   printf("#  SSqs     FPSqs            CPSqs            BPSqs             SOqs            FPOqs            CPOqs            BPOqs            CSLqs\n");

   double start = omp_get_wtime();

   
   int n = 300000;

   //float *v = malloc(n*sizeof(float));
   float *w = malloc(n*sizeof(float));

   double tSSqs  = 0;
   double tFPSqs = 0;
   double tCPSqs = 0;
   double tBPSqs = 0;
   double tSOqs  = 0;
   double tFPOqs = 0;
   double tCPOqs = 0;
   double tBPOqs = 0;
   double tCSLqs = 0;

   for(int r=1; r<=rounds; r++) {
      //fill(v, n);
      tSSqs  += etime(SSqs,  v, n, w, 0);
      tFPSqs += etime(FPSqs, v, n, w, 0);
      tCPSqs += etime(CPSqs, v, n, w, 0);
      tBPSqs += etime(BPSqs, v, n, w, 0);
      tSOqs  += etime(SOqs,  v, n, w, 0);
      tFPOqs += etime(FPOqs, v, n, w, 0);
      tCPOqs += etime(CPOqs, v, n, w, 0);
      tBPOqs += etime(BPOqs, v, n, w, 0);
      tCSLqs += etime(CSLqs, v, n, w, 0);
   }

   tSSqs  /= rounds;
   tFPSqs /= rounds;
   tCPSqs /= rounds;
   tBPSqs /= rounds;
   tSOqs  /= rounds;
   tFPOqs /= rounds;
   tCPOqs /= rounds;
   tBPOqs /= rounds;
   tCSLqs /= rounds;

   printf("%9.5f %9.5f %4.2f   %9.5f %4.2f   %9.5f %4.2f   %9.5f %4.2f   %9.5f %4.2f   %9.5f %4.2f   %9.5f %4.2f   %9.5f %4.2f\n", tSSqs, tFPSqs, tSSqs/tFPSqs, tCPSqs, tSSqs/tCPSqs, tBPSqs, tSSqs/tBPSqs, tSOqs, tSSqs/tSOqs, tFPOqs, tSSqs/tFPOqs, tCPOqs, tSSqs/tCPOqs, tBPOqs, tSSqs/tBPOqs, tCSLqs, tCSLqs/tBPOqs);

   //free(v);
   free(w);
   

   printf("\n# Total elapsed time for PERFORMANCE experiment: %5.2fs\n\n\n", omp_get_wtime()-start);
}

int main(void) {

  
        int n_hits = 300000;

        //Opening files

        //float dist[n_hits],x_bag[n_hits],y_bag[n_hits],z_bag[n_hits], x_pred, y_pred, z_pred;

        float x_pred, y_pred, z_pred;

        float *dist  = malloc(n_hits * sizeof(float));
        float *x_bag = malloc(n_hits * sizeof(float));
        float *y_bag = malloc(n_hits * sizeof(float));
        float *z_bag = malloc(n_hits * sizeof(float));


        int i = 0;
        FILE *fx, *fy,*fz;

        fx = fopen("input_datasets/bag_x.dat", "r");
        fy = fopen("input_datasets/bag_y.dat", "r");
        fz = fopen("input_datasets/bag_z.dat", "r");

        for(i=0; i < n_hits; i++){
                fscanf(fx, "%f", &x_bag[i]);
                fscanf(fy, "%f", &y_bag[i]);
                fscanf(fz, "%f", &z_bag[i]);
        }
        fclose(fx);
        fclose(fy);
        fclose(fz);


        x_pred = 13.0;
        y_pred = -28.0;
        z_pred = 12.78;

        #pragma omp parallel for
        for(i=0; i<n_hits;i++){
                double temp = pow((x_bag[i] - x_pred),2) + pow((y_bag[i] - y_pred),2) + pow((z_bag[i] - z_pred),2);
                dist[i] = sqrt(temp);
                //printf("%f \n",dist[i]);
        }

        for (i=299990; i < n_hits; i++){
                printf("x[%d] = %f\n", i, x_bag[i]);
                printf("y[%d] = %f\n", i, y_bag[i]);
                printf("z[%d] = %f\n", i, z_bag[i]);
                printf("dist[%d] = %f\n", i, dist[i]);
                printf("\n");
        }

        free(x_bag);
        free(y_bag);
        free(z_bag);
                        


	signal(SIGSEGV, error);
	srand(time(NULL));
	omp_set_nested(1);
	omp_set_num_threads(4);
	//find_cutoff(permutation, 10);
	performance(dist, 10);  
	return 0;
}

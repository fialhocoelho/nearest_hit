#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

#define N_HITS 300000
#define X_PRED 13.0f
#define Y_PRED -28.0f
#define Z_PRED 12.78f
#define SILENT false

#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

// Simple auxiliary inline functions

inline void swap(float v[], int i, int j) { float x=v[i]; v[i]=v[j]; v[j]=x; }
inline void index_swap(int v_index[], int i, int j) { int x=v_index[i]; v_index[i]=v_index[j]; v_index[j]=x; }
inline double lg(double x, double b) { return log(x)/log(b); }
//inline int best_at(int n) { return max(floor(lg(n,4)-4),1); }
//inline int best_cutoff(int n) { return max(n>>best_at(n),1); }

int best_at(int n){ 
	int aux = floor(lg(n,4)-4);
	return max(aux,1); 
}

int best_cutoff(int n) { 
	int aux = n>>best_at(n);
	return max(aux,1);
}


// Hoare partition algorithm (avoids worst case for sorted sequences)

int partition(float v[], int v_index[], int start, int end) {
   float pivot = v[(start + end)/2];
   int cut = end+1;
   start--;
   while( start < cut ) {
      do cut--;   while( v[cut]   > pivot );
      do start++; while( v[start] < pivot );
      if( start < cut ){ 
	      swap(v, start, cut);
	      index_swap(v_index, start, cut);
      }
   }
   return cut;
}


// Serial Standard Quicksort (may cause stack overflow)

void ssqs(float v[], int v_index[], int start, int end) {
   if( start >= end ) return;
   int cut = partition(v, v_index, start, end);
   ssqs(v, v_index, start, cut);
   ssqs(v, v_index,  cut+1, end);
}

void SSqs(float v[], int v_index[],int n) {
   ssqs(v, v_index, 0, n-1);
}



// PSqs - Parallel Standard Quicksort

void psqs(float v[], int v_index[], int start, int end, int cutoff) {
   if( end-start+1 <= cutoff )
      ssqs(v, v_index, start, end);
   else {
      int cut = partition(v, v_index, start, end);
      #pragma omp task
      psqs(v, v_index, start, cut, cutoff);
      #pragma omp task
      psqs(v, v_index, cut+1, end, cutoff);
   }
}

void PSqs(float v[], int v_index[], int n, int cutoff) {
   #pragma omp parallel
   #pragma omp single nowait
   psqs(v, v_index, 0, n-1, cutoff);
}

// SOqs  - Serial Optimized Quicksort

void soqs(float v[], int v_index[], int start, int end) {
   while( start < end ) {
      int cut = partition(v, v_index, start, end);
      if( cut-start+1 < end-cut ) {
         soqs(v,v_index, start, cut);
         start = cut + 1;
      }
      else {
         soqs(v, v_index, cut+1, end);
         end = cut;
      }
   }
}

void SOqs(float v[], int v_index[], int n) {
   soqs(v, v_index, 0, n-1);
}


// POqs - Parallel Optimized Quicksort

void poqs(float v[], int v_index[], int start, int end, int cutoff) {
   while( end-start+1 > cutoff ) {
      int cut = partition(v, v_index, start, end);
      if( cut-start+1 < end-cut ) {
         #pragma omp task 
         poqs(v, v_index, start, cut, cutoff);
         start = cut+1;
      }
      else {
         #pragma omp task
         poqs(v, v_index, cut+1, end, cutoff);
         end = cut;
      }
   }
   soqs(v, v_index, start, end); 
}

void POqs(float v[], int v_index[], int n, int cutoff) {
   #pragma omp parallel
   #pragma omp single nowait
   poqs(v, v_index, 0, n-1, cutoff);
}


// BPOqs - Best-grained Parallel Optimized Quicksort

void BPOqs(float v[], int v_index[], int n) {
   #pragma omp parallel
   #pragma omp single nowait
   poqs(v,v_index, 0, n-1, best_cutoff(n));
}



int main (int argc, char *argv[]){
	
        //omp_set_num_threads(atoi(argv[1]));
	
	//Opening files
	float *dist  = malloc(N_HITS * sizeof(float));
	float *x_bag = malloc(N_HITS * sizeof(float));
	float *y_bag = malloc(N_HITS * sizeof(float));
	float *z_bag = malloc(N_HITS * sizeof(float));
	

	int i = 0;
    	FILE *fx, *fy,*fz;
	int n_hits = N_HITS;
	double start_time, end_time, dist_time, sort_time;

    	fx = fopen("input_datasets/bag_x.dat", "r");
	fy = fopen("input_datasets/bag_y.dat", "r");
	fz = fopen("input_datasets/bag_z.dat", "r");

	for(i=0; i < N_HITS; i++){
		// if are used to ignore a fscan return value warning
		if(fscanf(fx, "%f", &x_bag[i]));
		if(fscanf(fy, "%f", &y_bag[i]));
		if(fscanf(fz, "%f", &z_bag[i]));
        }
        fclose(fx);
	fclose(fy);
	fclose(fz);


	start_time = omp_get_wtime();
	
	#pragma omp parallel for
	for(i=0; i<N_HITS;i++){
		double temp = pow((x_bag[i] - X_PRED),2) + pow((y_bag[i] - Y_PRED),2) + pow((z_bag[i] - Z_PRED),2);
		dist[i] = sqrt(temp);
		//printf("%f \n",dist[i]);
	}

	end_time = omp_get_wtime();
	dist_time = end_time - start_time;

        int n = 10;
        float *dist_for_serial;
        float *dist_for_parallel;

        float *v     = malloc(N_HITS * sizeof(float));
        int *v_index = malloc(N_HITS * sizeof(int));


        //memcpy(dist,   dist_for_serial, N_HITS * sizeof(float));
        //memcpy(dist_for_parallel, dist, N_HITS * sizeof(float));


	if (SILENT == true){
		for (i=299990; i < N_HITS; i++){
                	printf("x[%d] = %f\n", i, x_bag[i]);
                	printf("y[%d] = %f\n", i, y_bag[i]);
                	printf("z[%d] = %f\n", i, z_bag[i]);
                	printf("dist[%d] = %f\n", i, dist[i]);
			//printf("dist[%d] = %f\n", i, dist_for_serial[i]);
			//printf("dist[%d] = %f\n", i, dist_for_parallel[i]);
			printf("\n");
        	}
	}
	

	// Sort test
	// initialize indexes
	for (i=0; i< N_HITS; i++){v_index[i] = i;}	


	start_time = omp_get_wtime();
	
	//SSqs(dist,v_index, N_HITS);
        BPOqs(dist,v_index,N_HITS);

	end_time = omp_get_wtime();
	sort_time = end_time - start_time;
	
	//check sort results
	//for (i=0;i<n;i++){
	//	printf("%d->",v_index[i]);
        //      printf("%f ",dist[i]);
        //}
	//printf("\n\n%d\n",best_cutoff(N_HITS));
	
	free(x_bag);
        free(y_bag);
        free(z_bag);
	free(dist);
	//free(dist_for_serial);
	//free(dist_for_parallel);

	printf("dist_time: %lf \nsort_time: %lf\ntotal_time: %f\n", dist_time, sort_time, dist_time + sort_time);

	return 0;
}

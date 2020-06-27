#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <math.h>


int main (int argc, char *argv[]){
	int n=3, n_hits = 300000;

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
	

	return 0;
}

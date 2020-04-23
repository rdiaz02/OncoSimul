/*
	correlation functions at distance d=1 to d=max
	only binary landscapes
*/

#include "landscape.h"
#include "genotypes.h"
#include "summary_statistics.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>


int main(int argc, char *argv[]){

	int details=0;
	
	struct landscape h = ReadFile(argv[1], 0);
	
	int i,j;

	int l=h.nlocus;

	/*int v[10];
	  for(i=0;i<10;i++){
	  v[i]=2;
	  };
	  seed_ran1(-1);
	  for(n=1;n<=8;n++){
	  for(j=0;j<10;j++){
	  init_landscape(&h,10,v);
	  Kaufman_NK(&h,j,1,0);
	  int l=h.nlocus;
	  for(i=1;i<l;i++){
	  printf("gamma_dist_global %u %u %f\n",j,i,gamma_dist_global(h,i));
	  };
	  free_landscape(&h);
	  };
	  };*/
	  
	  
	h = ReadFile(argv[1], 0);
	  
	  
	/* work in log scale */

	for(i=0;i<h.ngenotypes;i++){
		h.fitness[i]=log(h.fitness[i]);
	}
	

	if(details==1){

		printf("gamma_den %f\n",gamma_den(h));

		for(i=1;i<=l;i++)
			printf("gamma_den_i (%u) %f\n",i,gamma_den_i(h,i));


		for(j=1;j<=l;j++)
			printf("gamma_den_j (%u) %f\n",j,gamma_den_j(h,j));
		
		
		printf("gamma_num_global %f\n",gamma_num_global(h));
		
		for(i=1;i<=l;i++)
			printf("gamma_num_i_global (%u) %f\n",i,gamma_num_i_global(h,i));
		
		
		for(j=1;j<=l;j++)
			printf("gamma_num_j_global (%u) %f\n",j,gamma_num_j_global(h,j));
		
		
		
		for(i=1;i<=l;i++)
			for(j=1;j<=l;j++)
				printf("gamma_num_ij_global (%u,%u) %f\n",i,j,gamma_num_ij_global(h,i,j));
		
		for(i=2;i<l;i++)
			printf("gamma_num_dist_global (d=%u) %f\n",i,gamma_num_dist_global(h,i));

		for(j=1;j<=l;j++)
			for(i=2;i<l;i++)
				printf("gamma_num_dist_j_global (d=%u,%u) %f\n",i,j,gamma_num_dist_j_global(h,i,j));


		printf("\n");
		printf("--------------\n");
		printf("\n");
	}
	
	printf("gamma_global %f\n",gamma_global(h));
	
	for(i=1;i<=l;i++){
		printf("gamma_i_global (%u) %f\n",i,gamma_i_global(h,i));
	}
	
	for(j=1;j<=l;j++){
		printf("gamma_j_global (%u) %f\n",j,gamma_j_global(h,j));
	}
	
	for(i=1;i<=l;i++)
		for(j=1;j<=l;j++)
			printf("gamma_ij_global (%u,%u) %f\n",i,j,gamma_ij_global(h,i,j));

	
	for(i=2;i<l;i++)
		printf("gamma_dist_global (d=%u) %f\n",i,gamma_dist_global(h,i));
	
	for(j=1;j<=l;j++)
		for(i=2;i<l;i++)
			printf("gamma_dist_j_global (d=%u,%u) %f\n",i,j,gamma_dist_j_global(h,i,j));


	/*  other stuff...
	    printf("gamma_global %f\n",gamma_global(h));
	    printf("gamma_i_global %f\n",gamma_i_global(h,1));
	    printf("gamma_i_global %f\n",gamma_i_global(h,2));
	    printf("gamma_j_global %f\n",gamma_j_global(h,1));
	    printf("gamma_j_global %f\n",gamma_j_global(h,2));
	    printf("gamma_ij_global %f\n",gamma_ij_global(h,1,2));
	    printf("gamma_dist_global %f\n",gamma_dist_global(h,1));
	    printf("gamma_dist_j_global %f\n",gamma_dist_j_global(h,1,1));
	    printf("gamma_dist_j_global %f\n",gamma_dist_j_global(h,1,2));  
	 */
	 
	return 0;
}

#include "generalized_chain.c"

int main(int argc, char *argv[]){

	struct landscape h;
	int i;
	
	int n=h.ngenotypes;

	double *ha0,*ha1;
	double *flux;
  	int *m;

	h=ReadFile(argv[1],0);
  	
  	flux=malloc(sizeof(double)*n*n);
  	m=malloc(sizeof(int)*n*n);
  	ha0=malloc(sizeof(double)*n);
  	ha1=malloc(sizeof(double)*n);
    

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
  
	/*
		work in log scale
	*/
	
	for(i=0;i<h.ngenotypes;i++)
		h.fitness[i]=log(h.fitness[i]);
  
  
	initialize_flux(h,flux);
	adjacency_matrix(h,m);
	assign_fluxes(h,m,flux);
	print_chain_pairs(h,flux,0.0001);
  
	/* 
		other stuff, correct but NOT WORKING

	get_node_constraints(h,m,ha0,ha1);  
	for(i=0;i<n;i++){
	  printf("H= %u %f %f\n",i,ha0[i],ha1[i]);
	}
	*/
  
	free(flux);
	free(m);
	free(ha0);
	free(ha1);
  
	return 1;
}


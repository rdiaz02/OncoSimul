/*
	Functions that count chains
	(genotypes with a single fitter neighbor that are 'chained')
	
	last update 18 nov 2013
	
*/

#include <stdlib.h>
#include <stdio.h>

#include "landscape.h"
#include "genotypes.h"
#include "summary_statistics.h"

/*
	For each genotype, compute the distribution
	of degree of outgoing steps
*/
int *compute_HistoBetterNeighbor( struct landscape * FL, float FitnessRatio )
{


	int q;

	int neighbors = 0;

	int *HistoBetterNeighbor;


	for( q=0; q<FL->nlocus; q++)
		neighbors += FL->alleles[q] - 1;


	/*
		Steps - Better Neighbors of each genotype
	*/
	HistoBetterNeighbor = (int *)calloc( (size_t) neighbors+1, (size_t) sizeof(int) );
	if( ! HistoBetterNeighbor  )fprintf(stderr, "compute_HistoBetterNeighbor: cannot allocate memory for HistoBetterNeighbor, bye\n"), exit(1);
	
	for(q=0;q < FL->ngenotypes;q++)
		HistoBetterNeighbor[  CountFitterNeighbors( FL, q, FitnessRatio, 1 ) ]++;

	return HistoBetterNeighbor;

}

/*
*/
void Connect_Fitness(struct FromTo * Ptr, struct landscape *fl, int g, float FitnessRatio){

	int l, a, g2;
	
	int *genotype = int2genotype( *fl, g , NULL);
	
	int tmp;
	
	for(l=0; l<fl->nlocus; l++){
	
		tmp = genotype[l];

		for(a=0; a< fl->alleles[l]; a++  ){
	
			if( a != tmp ){
				
				genotype[l]=a;
						
				g2 = genotype2int( *fl, genotype );
				
				if( fl->fitness[ g2 ] > fl->fitness[g]* FitnessRatio ){
		
					Ptr[g].To = g2;   /* connect the From->To */
					
					
					/*
						handle multiple origins
					*/
					if( Ptr[g2].From == NULL ){
						 Ptr[g2].From = (long *) malloc( (size_t) 2*sizeof(long)  );
						 if( ! Ptr[g2].From )fprintf(stderr, "Connect_Fitness: cannot allocate Ptr[%d].From, bye\n", g2), exit(3);
						 
						Ptr[g2].From[0] = 1;
						Ptr[g2].From[1] = g;
					}
					else{
						long x=Ptr[g2].From[0];
						
						Ptr[g2].From = (long *) realloc( Ptr[g2].From, (size_t) (x+2)*sizeof(long)  );
						if( ! Ptr[g2].From )fprintf(stderr, "Connect_Fitness: cannot reallocate Ptr[%d].From, bye\n", g2), exit(3);
						
						Ptr[g2].From[0] ++ ;
						Ptr[g2].From[ Ptr[g2].From[0] ] = g;
						
					}

				}
			
			}

			genotype[l] = tmp;
			
		}
	}

	free(genotype);

}


/*
	Extract recurssively the longest chain
*/
long GetDepth(  struct FromTo *Ptr, int g, struct landscape *fl ){

  long L=0;
  //long x,l, g_from;
  long x,l; 
	
  if( Ptr[g].From == NULL )
    return 0;

	
  //g_from = Ptr[g].From[1];
	
  for( x=0; x < Ptr[g].From[0]; x++ )
    {
		
      l =  GetDepth(  Ptr, Ptr[g].From[1+x] , fl) + 1;
		
      if( l>L )
	{
	  //g_from = Ptr[g].From[1+x];
	  L = (l>L)?l:L;
	}
    }
	

  return L;
}


/*
	Extract recurssively the sum of all chain steps in an independent chain
*/
long GetSteps( struct FromTo *Ptr, int g,  struct landscape *fl  ){

	long L=0;
	long x;
	
	if( Ptr[g].From == NULL )
		return 0;
	
	for( x=0; x < Ptr[g].From[0]; x++ )
		L +=  GetSteps(  Ptr, Ptr[g].From[1+x], fl ) + 1;


	return L;
}



/*
	Extract recursively the number of starting points
*/
long GetStartingPoints( struct FromTo *Ptr, int g, struct landscape *fl  ){

	long L=0;
	long x;
	
	if( Ptr[g].From == NULL )
		return 1;

	for( x=0; x < Ptr[g].From[0]; x++ )
		L +=  GetStartingPoints(  Ptr, Ptr[g].From[1+x], fl );

	return L;
}

/*
	Erase all connections
*/
void EraseChain( struct FromTo *Ptr, int g, struct landscape *fl  ){

	int x;

	if( Ptr[g].From == NULL ){
		Ptr[g].To = -1;        /* erase this node */
		return;
	}
	
	for( x=0; x < Ptr[g].From[0]; x++ )
		EraseChain(  Ptr, Ptr[g].From[1+x], fl );



	free( Ptr[g].From );
	Ptr[g].From = NULL;
	Ptr[g].To = -1; 

	return ;
}



/*
	fl is the fitness landscape

	When opt_sum is 0, take the longest section of all steps connected together
	                1, take the sum of all steps connected together
	
	IncreaseRatio is the minimum increase ratio (e.g. set to 1 for absolute increase)
*/
struct Chains Compute_Chains( struct landscape *fl , short opt_sum, float FitnessRatio){


	struct Chains my_chains; /* statistics of chains */
	
	long i;
	
	struct FromTo *Ptr;   /* all chains */
	
	long *Out1;           /* store all 1-way-out by their genotypes number */
	long nout1=0;
	
	struct FromTo *X;     /* this is used to ecxplore the current chain */
		
	
	Ptr = (struct FromTo *)malloc( (size_t)fl->ngenotypes*sizeof(struct FromTo) );
	Out1 = (long *)malloc( (size_t) fl->ngenotypes*sizeof(long) );
	
	if(!Ptr || !Out1 )fprintf(stderr, "ComputeChains: cannot allocate Ptr or Out1, bye\n"), exit(3);
	
	/*
		Init chains in the FromTo Ptr structure
	*/
	for(i=0;i<fl->ngenotypes; i++)
	{
		Ptr[i].To = -1;
		Ptr[i].From = NULL;
	}

	/*
		Built Chain Units
	*/
	for(i=0;i<fl->ngenotypes; i++){
	
		if(  CountFitterNeighbors(fl, i, FitnessRatio, 1 ) == 1 ){
			Out1[nout1++]=i;
			Connect_Fitness( Ptr, fl, i, FitnessRatio);
		}
	}
	
	/*
		Init the structure
	*/
	my_chains.nchains = 0;
	my_chains.steps	   = NULL;
	my_chains.depth	   = NULL;
	my_chains.origins = NULL;

	/*
		Conservative allocation
	*/
	if(nout1 != 0){

		my_chains.steps           = (long *) calloc( nout1, sizeof(long) );
		my_chains.depth           = (long *) calloc( nout1, sizeof(long) );
		my_chains.origins = (long *) calloc( nout1, sizeof(long) );
		if(!my_chains.steps || !my_chains.depth || !my_chains.origins)
			fprintf(stderr, "Compute_Chains: cannot allocate arrays, bye\n"), exit(3);
	}


	/*
		Count chains and their statistics
	*/
	
	for(i=0;i<nout1; i++)
	{
	
		if( Ptr[Out1[i]].From == NULL &&  Ptr[Out1[i]].To == -1 )                               /* it has been erased */
			continue;
		
		/*
			This genotype is part of a chain
		*/
		X = Ptr + Out1[i];               
		
		/*
			Follow it to the end-point
		*/
		while( X->To != -1 )
			X = Ptr + X->To; 
		
		/*
			Compute the number of spetsp (total or longest)
		*/
		my_chains.depth[ my_chains.nchains ]           = GetDepth(           Ptr , (long)(X - Ptr), fl );
		my_chains.steps[ my_chains.nchains ]           = GetSteps(           Ptr , (long)(X - Ptr), fl );
		my_chains.origins[ my_chains.nchains ] = GetStartingPoints(  Ptr , (long)(X - Ptr), fl );

		EraseChain( Ptr, (long)(X - Ptr), fl  );

		my_chains.nchains++;
		
	}
	
	/*
		Resize memory
	*/
	if(my_chains.nchains)
	{
		my_chains.steps           = (long *) realloc( (void *) my_chains.steps          ,(size_t) my_chains.nchains * sizeof(long) );
		my_chains.depth           = (long *) realloc( (void *) my_chains.depth          ,(size_t) my_chains.nchains * sizeof(long) );
		my_chains.origins = (long *) realloc((void *) my_chains.origins ,(size_t) my_chains.nchains * sizeof(long) );
		if(!my_chains.steps || !my_chains.depth || !my_chains.origins)
		fprintf(stderr, "Compute_Chains: cannot re-allocate arrays, bye\n"), exit(3);
	}
	
	/*
		A simple check 
	*/
	for(i=0;i<fl->ngenotypes; i++)
		if( Ptr[i].From )
			printf("ComputeChains_Fitness: this should not happen (chain still in memory)\n"),exit(1);



	free(Ptr);
	free(Out1);
	
	return my_chains;
}


void free_chain(  struct Chains mychain ){

	free( mychain.depth );
	free( mychain.steps );
	free( mychain.origins );

}

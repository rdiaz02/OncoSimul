#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "random.h"
#include "landscape.h"
#include "genotypes.h"
#include "summary_statistics.h"
#include "models.h"
#include "verbose.h"


/*
	Generate a landscape with complete epistasis --maximal dimansion
	though, only random noise
	equivalent in regards of fitness ranks to Kaufman_NK with K=N-1 or RMF, with f=1
	if sigma is <=0, use a Uniform[-0.5,+0.5]
*/

void HouseOfCards( struct landscape *l, float sigma){
  /* extern char verbose; */
	int i;
	
	if( sigma <= 0 )
		for( i=0; i<l->ngenotypes; i++)
		{
			l->fitness[i] = (float) uniform_dev( )-0.5;
		}
	else
		for( i=0; i<l->ngenotypes; i++)
		{
			l->fitness[i] = normal_dev(0, sigma );
		}


/*	if(opt_log)
		for( i=0; i<l->ngenotypes; i++)
			l->fitness[i] = exp( l->fitness[i] );
*/
		
	if(verbose)
		output_landscape(l);

}

/*
	It generates a landscape where each genotype has a multiplicative contribution
	-- additive in logscale --
	-- this corresponds to NK landscape with K=0 --

	l is the fitness landscape
*/
void Multiplicative( struct landscape *l , double s, double stdev, int  opt_single_fitness, float DimRet,  char opt_log )
{

  /* extern char verbose; */
	int *genotype=NULL;
	float *fa;  
	
	int g,n,
	a=0,
	t=0;
	
	do
		t += l->alleles[a++];
	while(a<l->nlocus);
	
	fa = (float *) malloc( (size_t) t*sizeof(float) );
	if(!fa)fprintf(stderr, "Multiplicative: cannot allocate the fitness of alleles, bye"), exit(3);

	if ( s==-1 )
	{
		float fit=0;
		
		if( opt_single_fitness == 1 )
			fit = (float) 1+uniform_dev( );       /* random [1,2] */
		
		for( a=0; a< t; a++){
			
			if (opt_single_fitness==0)
				fa[a] = (opt_log==0)?(float) 1+uniform_dev( ): uniform_dev( );  /* each mutation is unique */
			else
				fa[a] = (opt_log==0)? fit: fit-1;                       /* all mutation are identical */
		   		
			}
	}
	else if (s>-1 && stdev>0){
	
		for( a=0; a< t; a++){
			fa[a]= (opt_log==1)? normal_dev(s, stdev ) : normal_dev(1+s, stdev );
		}
			
	}
	else if ( s>-1 && stdev <=0 ){
	
		for( a=0; a< t; a++)
			fa[a]= (opt_log==0)?(float) 1.0+s: s;
	
	}
	else
	 	fprintf(stderr,"No model for your s along with its stdev\n"),exit(1);		 

	/*
		reset all 0 alleles at each locus to a contribution of 1
	*/
	t=0;a=0;	
	do
	{
		fa[t]=(opt_log)?0:1;                   /* no effect for the zero allele of each locus */
		t += l->alleles[a++];
	}
	while(a<l->nlocus);


	if( verbose )
	 	for( a=0; a< t; a++)
			printf("fa[%d]: %f\n", a, fa[a]);
	 

	for(g=0;g<l->ngenotypes; g++){

		int total_alleles=0;
		int dosage=0;


		genotype = int2genotype( *l , g, genotype );

		l->fitness[g]=(opt_log)?0:1;
		
	
		for(n=0; n< l->nlocus; n++){
		
			dosage += genotype[n];
			
			if(opt_log)
				l->fitness[g] += fa[ total_alleles+genotype[n] ];
			else
				l->fitness[g] *= fa[ total_alleles+genotype[n] ];
				
			total_alleles += l->alleles[n];
		}
		
		if(opt_log)
			l->fitness[g] += dosage*DimRet;
		else
			l->fitness[g] *= pow(1+DimRet, dosage);
		
		
	}
		
/*	if(opt_log)
		for( g=0; g<l->ngenotypes; g++)
			l->fitness[g] = exp( l->fitness[g] );
*/
	 	
}


/*
	It generates a landscape where there is an optimal number of mutated allele
	l is the fitness landscape
*/

void Optimum( struct landscape *l , double mean_p, double sigma_p, double mean_o, double sigma_o )
{

  /* extern char verbose; */
	int *genotype=NULL;
	float *fa;  
	
	int g,n,
	a=0,
	t=0;	


	do
		t += l->alleles[a++];
	while(a<l->nlocus);
	
	fa = (float *) malloc( (size_t) t*sizeof(float) );
	if(!fa)fprintf(stderr, "mean_pltiplicative: cannot allocate the fitness of alleles, bye"), exit(3);


	if ( mean_p==0 )
	{
		for( a=0; a< t; a++)
		{
			fa[a] = uniform_dev( );  /* each mean_ptation is unique */
		}
	}
	else if (mean_p>0 && sigma_p>0){
	
		for( a=0; a< t; a++)
			fa[a]= normal_dev(mean_p, sigma_p );	
			
	}
	else if ( mean_p>0 && sigma_p <=0 )
	{
	
		for( a=0; a< t; a++)
			fa[a]= mean_p;
	
	}
	else
	 	fprintf(stderr,"Optimum: No model for your sigma_p and mean_p\n"),exit(1);		 



	/*
		reset all 0 alleles at each locus to a contribution of 0
	*/
	t=0;a=0;	
	do
	{
		fa[t]=0;                   /* no effect for the zero allele of each locus */
		t += l->alleles[a++];
	}
	while(a<l->nlocus);

	if( verbose )
	 	for( a=0; a< t; a++)
			printf("prod[%d]: %f\n", a, fa[a]);
	 




	for(g=0;g<l->ngenotypes; g++){

		int total_alleles=0;
		double product=0;

		genotype = int2genotype( *l , g, genotype );
		
		for(n=0; n< l->nlocus; n++){
						
			product += fa[ total_alleles+genotype[n] ];
			total_alleles += l->alleles[n];
			
			l->fitness[g] = exp(  - pow( (mean_o - product),2)/(2*sigma_o) );
			
		}
		
		if( verbose ){
		 	print_intgenotype(g, l);
			printf("prod: %f, fit: %f\n", product, l->fitness[g] );
		}	 

		
	}
		
	free(fa);
}




#define PREVIOUS( n, L )  ((n)=((n)==0)?((L)-1):((n)-1))
/*
	Built a Kaufman NK landscapes
	K is the number of interacting loci ( K \in [0, N-1] )
	if random is 1, the K interacting loci are chosen at random ; 0 means neighboring sites
*/

void Kaufman_NK( struct landscape *l,  int K, short random ){


	/* extern char verbose; */

	int i,j, g;
	
	struct landscape *epistasis;
	int *sublandscape;
	
	long **random_loci=NULL;
	
	
	if( K > l->nlocus-1 || K<0){
		fprintf(stderr, "Kaufman_NK: K has to be 0 <= K < N\n");
		return;
	}
	

	epistasis = (struct landscape *)malloc( (size_t) l->nlocus*sizeof(struct landscape) );
	if(!epistasis)fprintf(stderr, "Kaufman_NK: cannot allocate epistasis, bye\n"), exit(3);
	
	
	sublandscape = (int *)malloc( (K+1)*sizeof( int )  );
	if(!sublandscape)fprintf(stderr, "Kaufman_NK: cannot allocate sublandscape, bye\n"), exit(3);
	
	
	if(random){
		random_loci = (long **)malloc( l->nlocus*sizeof(long *) );
		if(!random_loci)fprintf(stderr, "Kaufman_NK: cannot allocate random_loci\n"), exit(3);
	}
	
	
	/*
		This sets the epistasis for each locus
		assuming that interacting locis are neighbor
		half rounded on the left and half rounded-up
		on the right
	*/
	
	for( i=0; i< l->nlocus ; i++ ){
		
		/*
			copy sub-landscape
		*/
		if(! random ){
			
			int locus;
			
			locus = (i- K/2>=0)? i- K/2 : i- K/2 + l->nlocus;
			
			for( j=0;j<=K;j++ ){
				sublandscape[j] = l->alleles[locus];
				locus= ( locus == l->nlocus-1 )?0:locus+1;
			}
			
		}else{
		
			random_loci[i] = uniform_ArrayDiffInt(0, l->nlocus-1, K, i );

			sublandscape[0] = l->alleles[i];                                        /* The first one is the locus itself, the K others are randomnly chosen */

			for( j=0;j<K;j++ )
				sublandscape[j+1] = l->alleles[ random_loci[i][j] ];
		
			if(verbose){
				printf("locus %d; random: ", i);
				for( j=0; j<K; j++ )
					printf("%ld; ", random_loci[i][j]);
				printf("\n");
			}
		
		}
		
		
		
		/*
			init landscape
		*/
		
		init_landscape( epistasis+i, K+1, sublandscape);
		
		
		/*
			compute fitness
		*/
		
		HouseOfCards( epistasis+i ,-1 );
		
		if(verbose){
			printf(">> LOCUS %d\n", i);
			print_landscape( epistasis+i );
		}
		
	}
	
	
	/*
		This combines the epistasis of each locus into
		a whole fitness landscape
	*/

	for(g=0;g<l->ngenotypes;g++){
	
		int *genotype = int2genotype( *l, g , NULL);
//		l->fitness[g] = 0;
		l->fitness[g] = 1;
		for( i=0; i < l->nlocus; i++ ){
				
			/*
				copy sub-genotype
			*/

			if(! random ){
			
				int locus;
			
				locus = (i- K/2>=0)? i- K/2 : i- K/2 + l->nlocus;
			
				for( j=0;j<=K;j++ ){
					sublandscape[j] = genotype[locus];
					locus= ( locus == l->nlocus-1 )?0:locus+1;
				}
			
			}else{
		
				sublandscape[0] = genotype[i];                                        /* The first one is the locus itself, the K others are randomnly chosen */

				for( j=0;j<K;j++ )
					sublandscape[j+1] = genotype[ random_loci[i][j] ];
		
			}
			
			/*
				Extract and add corresponding fitness
			*/
			l->fitness[g] += 0.5+epistasis[i].fitness[ genotype2int( epistasis[i] , sublandscape ) ];
		}
		
		
		free(genotype);
	}


	
	if(random){
		for(i=0; i<l->nlocus; i++) free(random_loci[i]);
		free(random_loci);
	}

	free(sublandscape);

	for(i=0; i<l->nlocus; i++) free_landscape( epistasis+i );
	free(epistasis);
}



/*
	Each locus interacts with its neighbors
*/
void Ising( struct landscape *l, float mu_c, float sigma_c, short opt_circ )
{
  /* extern char verbose; */
	int i,g;
	
	float *costs;
	
	
	costs = (float *) malloc( (size_t) l->ngenotypes*sizeof(float) );
	if(!costs)fprintf(stderr, "Ising: cannot allocate the fitness of alleles, bye"), exit(3);
    
	if (mu_c==0){
		
		for( i=0; i< l->nlocus ; i++ )
			costs[i] =  uniform_dev();
	
	}
	else if (sigma_c<=0 && mu_c>0)
	{
	    for( i=0; i< l->nlocus ; i++ )
		costs[i] = mu_c;
	}
	else
		if (sigma_c>0 && mu_c>0)
		{
			for( i=0; i< l->nlocus ; i++ )
			costs[i] = normal_dev(mu_c, sigma_c );
		}
		
		
	if(verbose)
	{
		for( i=0; i< l->nlocus ; i++ )
			fprintf(stderr, "c[%d]= %f\n",i, costs[i]);
	}
		
	/*
		pour chaque genotype on regarde la compatibilite avec son voisin 'de droite'
	 */		
	for (g=0;g<l->ngenotypes;g++)
	{
		int *genotype = int2genotype( *l, g , NULL);
		
		double f=0;
		
		for( i=0; i < l->nlocus - 1; i++ )
		{
			if( genotype[i] != genotype[i+1] )
				f -= costs[i];
				
		}
				
		if( opt_circ )
			if (genotype[0]!=genotype[l->nlocus -1])
				f -= costs[l->nlocus-1];
		
			
		l->fitness[g]=f;
	}
	
	free(costs);
}






/*
	For now the simplest one
*/
void EggBox( struct landscape *l, float mu_e, float sigma_e ){


	int g;  /* iterate over all genotype */
	int i;  /* iterate over all locus */
	
	float effect;

	for(g=0;g<l->ngenotypes;g++){
	
		int *genotype = int2genotype( *l, g , NULL);
		int n1=0;
		
		for(i=0;i<l->nlocus;i++)
			n1 += genotype[i];
		
		if(sigma_e <= 0)
			effect = mu_e;
		else
			effect =  normal_dev(mu_e, sigma_e );
		
		l->fitness[g] = (n1%2)?effect/2:-effect/2;

	}

	
}





struct model_opt init_model( void  ){

	struct model_opt myoptions;
	
	myoptions.fix = 0;
	
	myoptions.kaufman_K = -1;
	myoptions.kaufman_rand = 0;
	
	myoptions.mu_s = -2;
	myoptions.sigma_s = -1;
	
	myoptions.sigma_hoc = -1;
	
	myoptions.mu_eggbox = -1;
	myoptions.sigma_eggbox = -1;
	
	myoptions.mu_ising = -1;
	myoptions.sigma_ising = -1;
	myoptions.circular_ising = 0;
	
	myoptions.DimRet = 0;
	
	myoptions.mu_optimum = -1;
	myoptions.sigma_optimum = -1;
	
	myoptions.mu_prod = -1;
	myoptions.sigma_prod = -1;

	return myoptions;

}


/*
	Values for all models are summed into a final fitness
	-- so log scale is probably the right scale.
*/
void MixedModel( struct landscape *fl, struct model_opt myoptions, int opt_log ){


	struct landscape fl_tmp;
	/* extern char verbose; */

	
	int g;
	




	/*
		The fixed amount
	*/
	for(g=0; g<fl->ngenotypes; g++)
		fl->fitness[ g ] = myoptions.fix;

	
		
	/*
		Multiplicative
	*/
	if( myoptions.mu_s > -1 )
	{
	
		if(verbose)
			printf("Mult mu_s: %f sigma_s: %f opt_log: %d\n", myoptions.mu_s, myoptions.sigma_s, opt_log);

		init_landscape( &fl_tmp, fl->nlocus, fl->alleles);
		Multiplicative( &fl_tmp, myoptions.mu_s, myoptions.sigma_s, 0, myoptions.DimRet, opt_log);

		for(g=0; g<fl->ngenotypes; g++)
			fl->fitness[ g ] += fl_tmp.fitness[ g ];

		free_landscape(&fl_tmp);
	}

	/*
		House of Cards
	*/
	if( myoptions.sigma_hoc >= 0 )
	{
	
		if(verbose)
			printf("HoC sigma_hoc %f\n", myoptions.sigma_hoc );
	
		init_landscape( &fl_tmp, fl->nlocus, fl->alleles);
		HouseOfCards(   &fl_tmp, myoptions.sigma_hoc);

		for(g=0; g<fl->ngenotypes; g++)
			fl->fitness[ g ] += fl_tmp.fitness[ g ];

		free_landscape(&fl_tmp);
	}
	
	/*
		Kaufmann NK
	*/
	if( myoptions.kaufman_K >=1 && myoptions.kaufman_K < fl->nlocus )
	{
	
		if(verbose)
			printf("NK\n");
	
		init_landscape( &fl_tmp, fl->nlocus, fl->alleles);
		Kaufman_NK(     &fl_tmp,  myoptions.kaufman_K, myoptions.kaufman_rand );

		for(g=0; g< fl->ngenotypes; g++)
			fl->fitness[ g ] += fl_tmp.fitness[ g ];

		free_landscape(&fl_tmp);

	}
	

	
	/*
		Eggbox
	*/
	if( myoptions.mu_eggbox > -1 )
	{
	
		if(verbose)
			printf("EggBox\n");
	
		init_landscape( &fl_tmp, fl->nlocus, fl->alleles);
		EggBox( &fl_tmp, myoptions.mu_eggbox, myoptions.sigma_eggbox );

		for(g=0; g< fl->ngenotypes; g++)
			fl->fitness[ g ] += fl_tmp.fitness[ g ];


		free_landscape(&fl_tmp);
	}


	/*
		Ising
	*/
	if( myoptions.mu_ising > 0 )
	{
	
		if(verbose)
			printf("Ising\n");
	
		init_landscape( &fl_tmp, fl->nlocus, fl->alleles);
		Ising( &fl_tmp, myoptions.mu_ising, myoptions.sigma_ising, myoptions.circular_ising );

		for(g=0; g<fl->ngenotypes; g++)
			fl->fitness[ g ] += fl_tmp.fitness[ g ];

		free_landscape(&fl_tmp);
	}
	

	
	/*
		Optimal product
	*/
	if( myoptions.mu_optimum > 0 && myoptions.mu_prod > 0 )
	{
	
		if(verbose)
			printf("Optimum\n");
	
		init_landscape( &fl_tmp, fl->nlocus, fl->alleles);
		Optimum( &fl_tmp , myoptions.mu_prod, myoptions.sigma_prod, myoptions.mu_optimum, myoptions.sigma_optimum );

		for(g=0; g<fl->ngenotypes; g++)
			fl->fitness[ g ] += fl_tmp.fitness[ g ];

		free_landscape(&fl_tmp);
	}
	
	if( opt_log )
	{
		fl->log_scale=1;
	}

		
}










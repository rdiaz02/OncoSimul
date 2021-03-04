/*
	From a starting genotype g, this functions
	count fitter genotypes after s mutations
	
	last update : nov, 18, 2013

*/


#include "landscape.h"
#include "genotypes.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <string.h>

/*
	Instantiate a new empty list
*/
struct list new_empty_list(){

	struct list l;
	l.n=0;
	l.genotypes=NULL;
	
	return l;
}

/*
	Free the list
*/
void free_list( struct list *l  ){

	if(l->genotypes)
		free(l->genotypes);
	l->genotypes = NULL;
	l->n=0;
}


/*
	if inc>0, it increases size list
	if <0, it shrinks the list
	if size is 0, free the list
*/
void resize_list( int inc, struct list *l ){

	if( l->n + inc == 0 ){
		free_list( l );
	}
	else{

		l->genotypes = (long *)realloc(  (void *)l->genotypes, (size_t) (l->n + inc)*sizeof(long) );
		if( ! l->genotypes )fprintf(stderr, "resize_list: cannot re-allocate l->gentoypes  from %ld to %ld, bye\n",  l->n , l->n + inc), exit(3);
	
		l->n += inc;
	}

}



/*
	Exchange both list
*/
void swap_list( struct list *l1, struct list *l2 ){

	struct list tmp;
	
	tmp = *l1;
	*l1 = *l2;
	*l2 = tmp;
}

/*
	return 1, 0 , -1 that corresponds to ( f2 larger ),  (both eq), ( f1 larger )
	the operator are modulated by the Fitness_Increase factor (set to 1 for absolute inc/dec/neutral)
*/
int compare_fitness( float f1, float f2, float Fitness_Increase ){

	if(f1 == DEFAULT_FITNESS || f2 == DEFAULT_FITNESS)
		fprintf(stderr, "compare_fitness: cannot compare fitness of an undefined genotype, bye\n"), exit(4);

	if( f2>f1 && f2>Fitness_Increase*f1)   /* f2 is really bigger than f1 */
		return 1;

	if( f1>f2 && f1>Fitness_Increase*f2 )   /* f1 is bigger */
		return -1;

	return 0;
}

/*
	return the number of neighbor with a all/better/worst fitness
	fl is the fitness landscape
	g is the genotype from which we look for neighbors	

	return the number of all/better/worst genotypes that is not already visited (info given by visited_genotypes)
	if visited_genotypes == NULL, ignore the 'visited' status

	the genotypes themselves are stored in the list li
	
	if opt_fit == "wnf" retrieve all genotypes (fitter, neutral, worst)
	if opt_fit =="f" retrieve fitter genotypes 
	if opt_fit == "w" retrieve worst genotypes
	if opt_fit == "n" retrieve neutral genotypes
	 
	Count as Fitter / Less Fit ONLY if the ratio is strictly larger than 'Fitness_Increase' (e.g. use 1 to get all Fitter genotypes)
	
*/

int RetrieveNeighbors( struct landscape *fl, int g, char *visited_genotypes, struct list *li, char *opt_fit, float Fitness_Increase ){


	int x=0,          /* number of neighboring genotypes */
	    nneighbours=0;    /* number of new neighboring genotypes - the ones that are stored */

	int *genotype = int2genotype( *fl, g , NULL);  /* the array form of the genotype --usefull for mutations */
	
	int l,            /* locus */
	    a,            /* allele */
	    g2,           /* neighbor */
	    tmp;          /* used to built neighbors (g2) */
	
	
	/*
		Quickly count how many are better/worst disregarding the visited status
	*/
	for( l=0; l < fl->nlocus; l++ ){
	
		tmp = genotype[l];
		
		for(a=0; a< fl->alleles[l]; a++){
		
			if( a != tmp ){
						
				genotype[l] = a;   /* mutate the locus into a new allele --> genotype is a neighbor */
				g2 = genotype2int( *fl, genotype);
				
				if (fl->fitness[g2] == DEFAULT_FITNESS)
					continue;
				
				if (strchr(opt_fit,'f')!=NULL)
					if( compare_fitness(fl->fitness[g] , fl->fitness[g2] , Fitness_Increase ) == 1 )
						x++;
				if (strchr(opt_fit,'w')!=NULL)
					if( compare_fitness(fl->fitness[g] , fl->fitness[g2] , 1.0/Fitness_Increase ) == -1 )
						x++;		
				if (strchr(opt_fit,'n')!=NULL)
					if( compare_fitness(fl->fitness[g] , fl->fitness[g2] , Fitness_Increase ) == 0 )
						x++;		
						
			}
		}
		genotype[l] = tmp;
	}
		

/*	printf("first scan is %d (opt_LessFit= %d)\n", x, (int)opt_LessFit); */
	
	resize_list( x, li );       /* resize list to its maximum value */
	
	
	/*
		now really check how many different genotypes there are and store them in list
	*/
	for( l=0; l < fl->nlocus; l++ ){
	
		tmp = genotype[l];
		
		for(a=0; a< fl->alleles[l]; a++){
		
			if( a != tmp ){

				genotype[l] = a;   /* mutate the locus into a new allele --> genotype is a neighbor */
				g2 = genotype2int( *fl, genotype);
				
				if (fl->fitness[g2] == DEFAULT_FITNESS)
					continue;
				
				if (strchr(opt_fit,'f')!=NULL)
					
					if( ( compare_fitness(fl->fitness[g] , fl->fitness[g2] , Fitness_Increase ) == 1) && (!visited_genotypes || ! visited_genotypes[ g2 ] ) ){
						li->genotypes[ li->n - x + nneighbours ] = g2; // li->n - x because sometimes li not empty! 
						if(visited_genotypes)
							visited_genotypes[ g2 ] = 1;
						nneighbours++;
					}
				
					
				if (strchr(opt_fit,'w')!=NULL)
					if( ( compare_fitness(fl->fitness[g] , fl->fitness[g2] , Fitness_Increase ) == -1) && (!visited_genotypes || ! visited_genotypes[ g2 ] ) ){
						li->genotypes[ li->n - x + nneighbours ] = g2;
						if(visited_genotypes)
							visited_genotypes[ g2 ] = 1;
						nneighbours++;
					}
	
				if (strchr(opt_fit,'n')!=NULL)
					if( ( compare_fitness(fl->fitness[g] , fl->fitness[g2] , Fitness_Increase ) == 0) && (!visited_genotypes || ! visited_genotypes[ g2 ] ) ){
						li->genotypes[ li->n - x + nneighbours ] = g2;
						if(visited_genotypes)
							visited_genotypes[ g2 ] = 1;
						nneighbours++;
					}



			}
		}
		genotype[l] = tmp;
	}



	if( x > nneighbours )                         /* if the list was too large resize it to the right size */
		resize_list( nneighbours-x, li );

	
	free(genotype);

	return nneighbours; 
}


/*
	g is the starting genotype
	if opt_cumul=1, cumulate the number of visisted genotypes, otherwise report their current number (after n steps)
	MaxSteps, stop the counting after MaxSteps. If set 0, then stop when the paths have all reached a peak.
	it returns the number of genotypes reached after s steps (from 1 to MaxSteps)	
*/
struct list count_genotypes( int g, struct landscape *fl, char opt_cumul, char opt_LessFit, int MaxSteps, float FitnessRatio ){

	char *visited_genotypes;  /* this store for all geneotytpes if it has been visited -- to handle redundancy -- */
	
	struct list l,            /* the main list */
	            l2,           /* the genotypes visited this round */
		    lcumul;       /* only used in case of cumulated genotypes (if you sums all new genotypes discovered so far */
		    
	struct list results;      /* a struture to store the results */


	int s=0,                  /* number of steps */
	    i,
	    new=0,
	    max;                 /* the maximum number of results to be stored. Equals MaxSteps if it is different from 0. Otherwise 10 by default and then increased by power of 2 */
	   
	int peak=0;              /* number of peaks at each step, s */
	

	/*
		init empty lists
	*/
	l       = new_empty_list();
	l2      = new_empty_list();
	results = new_empty_list();
	
	if(opt_cumul)
		lcumul = new_empty_list();

	
	
	max = ( MaxSteps )? MaxSteps : 10;
	resize_list( max, &results );
	
	/*
		seed list with g
	*/
	resize_list( 1, &l );
	l.genotypes[0]=g;
	
	if(opt_cumul){
		resize_list( 1, &lcumul );
		lcumul.genotypes[0]=g;
	}
	
	
	/*
		init visited_genotypes to zeros
	*/
	visited_genotypes = (char *)calloc( (size_t) fl->ngenotypes, (size_t) sizeof(char) );
	if(!visited_genotypes)fprintf(stderr, "count_genotypes: cannot c-allocate visited_genotypes, bye\n"), exit(3);
	
		
	while ( (! MaxSteps || s<MaxSteps) && l.n - peak > 0){
	
	
		/* for all genotype in the list */
			/* retrieve the fitter genotypes & add them to the new_list */
		
		peak=0;                    /* before evaluation at this step, no peak is reached */
		
		for( i=0; i<l.n; i++ ){
		
			int X;
			char opt_Fit[3];
			//int opt_Fit=(opt_LessFit==1)?-1:1;
			if (opt_LessFit==1) strcpy(opt_Fit,"f"); else strcpy(opt_Fit,"w");
		
			X = RetrieveNeighbors( fl, l.genotypes[i], visited_genotypes, &l2, opt_Fit, FitnessRatio );    /* this counts the number of 'new' genotypes with s mutations from g */
			
			
			if( X == 0  &&  visited_genotypes[ l.genotypes[i] ] == 0 ){                   /* this would be a case of a peak that is not already seen at this step */
				
				resize_list( 1 , &l2 );
				l2.genotypes[ l2.n - 1 ] = l.genotypes[i];
				visited_genotypes[ l.genotypes[i] ] = 1;
				
				peak++;
			}
			
			
		}

		
		/* if opt_cumul */
			/* merge list and list_2 while resetting unvisted_genotypes */
		/* else */
			/* replace list by list_2 while resetting and using unvisted_genotypes */


		for( i=0; i<l2.n; i++ )
			visited_genotypes[ l2.genotypes[i] ] = 0;
			
		if( opt_cumul ){
		
			
			new=0;                                                   /* true number of unvisited genotype so far */
		
			for( i=0; i<lcumul.n; i++ )
				visited_genotypes[ lcumul.genotypes[i] ] = 1;
			
			resize_list( l2.n , &lcumul );                           /* guessing that there will be n2 new genotype (actually, it is a maximum) */
			
			for( i=0; i<l2.n; i++ ){
			
				if( visited_genotypes[ l2.genotypes[i] ] == 0 ){
				
					lcumul.genotypes[ lcumul.n -l2.n + new ] = l2.genotypes[i];
					new ++;

				}
			}
			
			if( new<l2.n )
				resize_list( new - l2.n , &lcumul );

			for( i=0; i<lcumul.n-new ; i++ )
				 visited_genotypes[ lcumul.genotypes[i] ] = 0;
			
		}
		
		
		if(opt_cumul)
			results.genotypes[s] = lcumul.n;
		else
			results.genotypes[s] = l2.n;
		

		/*
			free l, copy l2 into l, and reset l2
		*/
						
		swap_list( &l, &l2 );
		free_list( &l2 );

		if( s >= max-1 ){
			
			resize_list( max, &results );
			max *= 2;
		}

		s++;		

	}
	
	
	
	if( max > s ){
		resize_list( s-max, &results );
	}
	
	free_list( &l );
	if(opt_cumul)free_list( &lcumul );
	
	free(visited_genotypes);


	return results;
}


/*
	return the number of genotpyes with a better fitness
	g is the genotype -- in int --
	IncreaseRatio is the minimum increase ratio (e.g. set to 1 for absolute increase)
	if opt_strict is set to 1 use the ">" operator, otherwise, use ">="
*/
int CountFitterNeighbors( struct landscape *fl, int g, float IncreaseRatio, int opt_strict ){


	int *genotype = int2genotype( *fl, g , NULL);
	int tmp;
	
	int nbetter=0;
	
	int l,     /* current locus */
	    a,     /* current allele */
	    g2;    /* neighbor */

	int cmp;   /* comparison 1, 0, -1 */


	for( l=0; l < fl->nlocus; l++ ){               /* for all loci */
		
		tmp = genotype[l];
		
		for(a=0; a< fl->alleles[l]; a++){     /* for all alleles */
			
			if( a != tmp ){

				genotype[l] = a;
				
				g2 = genotype2int( *fl, genotype);
		
				if (fl->fitness[g2] == DEFAULT_FITNESS)
					continue;
		
				cmp = compare_fitness( fl->fitness[g], fl->fitness[g2], IncreaseRatio);

				
				if( (cmp == 1 && opt_strict) ||  (cmp>=0 && !opt_strict) )
					nbetter++;
			}
			
		}
		
		genotype[l] = tmp;
	}


	free(genotype);

	return nbetter;
}


int CountFitterGenotypes( struct landscape *fl, int g, float IncreaseRatio, int opt_strict )
{
	
  //int tmp, nbetter=0;
  int nbetter = 0;
	int g2;     /* other genotype */
	int cmp;   /* comparison 1, 0, -1 */
		
	for(g2 = 0; g2 < fl->ngenotypes; g2++)  /* for all genotypes */
	{
		
		if( g != g2 )
		{
			
			if (fl->fitness[g2] == DEFAULT_FITNESS)
			{
				continue;
			}
			
			cmp = compare_fitness( fl->fitness[g], fl->fitness[g2], IncreaseRatio);
			if( (cmp == 1 && opt_strict) ||  (cmp>=0 && !opt_strict) )
			{
				nbetter++;
			}
		}
	}
	
	return (nbetter);
}

/*
	How many neighbor are defined (where the fitness value is not DEFAULT_FITNESS)
*/
int CountDefinedNeighbor( struct landscape *fl, int g){


	int *genotype = int2genotype( *fl, g , NULL);
	int tmp;
	
	int ndefined=0;
	
	int l, a;


	for( l=0; l < fl->nlocus; l++ ){
		
		tmp = genotype[l];
		
		for(a=0; a< fl->alleles[l]; a++){
			
			if( a != tmp ){

				genotype[l] = a;
				
				if(  fl->fitness[ genotype2int( *fl, genotype)] != DEFAULT_FITNESS )
					ndefined++;
			}
			
		}
		
		genotype[l] = tmp;
	}


	free(genotype);

	return ndefined;
}

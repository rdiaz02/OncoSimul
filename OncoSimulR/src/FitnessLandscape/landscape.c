/*

	Generate a Fitness Landscape on which each genotype (that is encoded by a long int)
	has an associated fitness value

*/

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

void print_genotype( int *genotype, int nlocus ){
	
	int i;

	for(i=0;i<nlocus-1; i++)
		printf("%d ", genotype[i] );
	
	printf("%d", genotype[nlocus-1] );
	

}
void print_intgenotype(int g, struct landscape *h){

	int *genotype=int2genotype( *h, g , NULL);
	print_genotype( genotype , h->nlocus );
	free( genotype );

}

int *int2genotype( struct landscape h, int x , int *genotype)
{
	int l=h.nlocus-1;
	
	if(genotype == NULL){
		genotype = (int *)malloc( (size_t) h.nlocus*sizeof(int) );
		if( ! genotype )fprintf(stderr, "int2genotype: cannot allocate the genotype, bye\n"), exit(3);
	}

	while( l >= 0 ){
		genotype[l] = x%h.alleles[l];
		
		x = x/h.alleles[l];
		l--;
	}

	return genotype;
}

/*
	built the index number that corresponds to the genotype
*/
int genotype2int( struct landscape h, int *genotype)
{


	int l=h.nlocus-1;
	int m=1;
	int x=0;
		
	while( l >= 0 ){
		
		x += genotype[l]*m;
		m*=h.alleles[l];
		l--;
	}

	return (x);
}

int isvalueinarray(int val, int *arr, int size)
{
	int i;
	for (i=0; i < size; i++)
	{
		if (arr[i] == val)
		{
			return (1);
		}
	}
	return (0);
}



void init_landscape( struct landscape *l, int nlocus, int *alleles)
{

	int i;

	l->nlocus = nlocus;
	l->neighbors = 0;
	l->nalleles = 0;
	l->alleles = (int *)malloc( (size_t) nlocus*sizeof(int) );
	l->nalleles = 0;
	
	if(! l->alleles )fprintf(stderr, "init_landscape: cannot allocate alleles, bye\n"), exit(3); 
	
	for(i=0;i<nlocus; i++)
	{
		l->alleles[i] = alleles[i];
		l->nalleles += alleles[i];
	}
	
	l->ngenotypes=1;
	for(i=0;i<l->nlocus;i++)
		l->ngenotypes *= l->alleles[i];
	
	l->fitness = (float *) malloc( (size_t) l->ngenotypes*sizeof(float) );
	if( !l )fprintf(stderr, "cannot allocate fitness\n" ), exit(3);
	
	for (i=0;i<l->ngenotypes;i++)
		l->fitness[i] = DEFAULT_FITNESS; /* it should be no fitness but be carefull... */
		
	for(i=0;i<l->nlocus;i++)
		l->neighbors += (l->alleles[i]-1);
		
	l->log_scale = 0;
}



void print_landscape( struct landscape *l )
{

	int g;
	int *geno=NULL;

	printf("landscape (%d loci; %d genotypes): ", l->nlocus, l->ngenotypes);
	print_genotype( l->alleles, l->nlocus );
	printf("\n");

	for(g=0; g<l->ngenotypes; g++ ){
	
		geno = int2genotype( *l, g , geno );
		printf("genotype: %d '", g);
		print_genotype( geno, l->nlocus );
		printf("'  %f\n",  l->fitness[g] );

	}
	printf("neighbors:%d\n",l->neighbors);
	printf("max=%f min=%f\n",l->minf,l->maxf);

	free(geno);
}

void output_landscape( struct landscape *l )
{

	int g;
	int *geno=NULL;
	
	print_genotype( l->alleles, l->nlocus );
	printf("\n");

	for(g=0; g<l->ngenotypes; g++ ){
	
		geno = int2genotype( *l, g , geno );
		print_genotype( geno, l->nlocus );
		printf(" %f\n",  l->fitness[g] );

	}

	free(geno);
}



void free_landscape( struct landscape *l )
{
	free(l->alleles);
	free(l->fitness);
}


void log_landscape( struct landscape *fl )
{

	int g;
	
	/*
		Basic Check
	*/
	for(g=0;g< fl->ngenotypes; g++ )
	{
		if(fl->fitness[g] <= 0 && fl->fitness[g]!=DEFAULT_FITNESS)
		{
			printf("<BR><BR>log_landscape: the FL contains <=0 (fitness=%f) fitness values, cannot use log(fitness)\n<BR>",fl->fitness[g]);
			print_intgenotype(g, fl);
			exit (1);
		}
	}
	
	/*
		Change into logscale
	*/
	for(g=0;g< fl->ngenotypes; g++ )
	{
		if (fl->fitness[g]!=DEFAULT_FITNESS)
		{
			fl->fitness[g] = log( fl->fitness[g] );
		}
	}
	
	fl->log_scale++;
	
	return;
}

void exp_landscape( struct landscape *fl )
{

	int g;
	
	/*
		Change into exp scale
	*/
	for(g=0;g< fl->ngenotypes; g++ )
	{
		 fl->fitness[g] = exp( fl->fitness[g] );
	}

	fl->log_scale--;

	
	return;
}

/*
	w is the weight of the first one
*/
struct landscape Merge( struct landscape *fl1, struct landscape *fl2, float w )
{

	struct landscape res;
	int i;
	
	if (w<0 || w>1)
	{
		fprintf(stderr,"weight should be between 0 and 1\n"),exit(1);
	}
		
	if( fl1->ngenotypes != fl2->ngenotypes )
	{
		fprintf(stderr,"cannot merge FL of different sizes\n"),exit(1);
	}


	init_landscape( &res, fl1->nlocus, fl1->alleles);
	
	/*
		Combine using weights in logscale
	*/
	for (i=0;i<res.ngenotypes;i++)
	{
		res.fitness[i] =  w*fl1->fitness[i] + (1-w)*fl2->fitness[i];
	}

	return (res);
}


void remove_negative( struct landscape *fl, char opt_cut )
{

	float min=fl->fitness[0];
	int i;

	if(opt_cut == 0)
	{
	
		for(i=1;i<fl->ngenotypes;i++)
		{
			if(fl->fitness[i]<min)
			{
				min=fl->fitness[i]<min;
			}
		}

		if(min<0)
		{
			for(i=0;i<fl->ngenotypes;i++)
			{
				fl->fitness[i] -= min;
			}
		}
		
	}
	else
	{
		for (i=0;i<fl->ngenotypes;i++)
		{
			if( fl->fitness[i] <= 0)
			{
				fl->fitness[i] = 0;
			}
		}
	}
	
}


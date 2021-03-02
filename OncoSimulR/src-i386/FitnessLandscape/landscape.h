/*
	This file lists all functions that deal with FL in general
	
	for more details on FL information extraction see
		genotypes.h
	
	for details about statistics on FL see
		summary_statistics.h
*/


#ifndef _LANDSCAPE_H_
#define _LANDSCAPE_H_

#define DEFAULT_FITNESS -99999999999

/* extern char verbose; */
/* if set to non-0 , become verbose */


/*
	Structure of a landscape
*/
struct landscape {

	int nlocus;      /* number of loci */
	int *alleles;    /* array with the number of alleles on each locus */
	int nalleles;	 /* number of alleles */
	
	int ngenotypes;  /* number of possible genotypes */
	float *fitness;  /* the fitness of each gentoype */
	
	char log_scale;      /* log-fitness indicator */ /* ::SHOULDN'T THIS BE AN INT:: */

	
	// :: ADD NUMBER OF OPTIMA AS MEMBER ? :: //
	
	int neighbors;   /* number of neighbors */
	float minf;      /* lowest fitness */
	float maxf;      /* highest fitness */
};



/*
	input.c
*/
struct landscape ReadFile( char *filename, int opt_zero );   /* if opt_zero is set to 1, fill missing fitness with 0 */
int * char_arg2array_int(  char *char_arg , int *nlocus );   /* change an char arg into an int array */
int * GetOrderFromFile( char *filename );                    /* associate each genotype to its rank in input file */

/*
	landscape.c
*/


/*
	print out a genotype
	int *genotype is the array form of a genotype
	nclous is the number of locus
*/
void print_genotype( int *genotype, int nlocus );
void print_intgenotype(int g, struct landscape *h);


/*
	Change a genotype into an array from an int or vice-versa
	
	if int *genotype == NULL when calling int2genotype,
		--> get the memory and return the array.
	Otherwise,
		--> simply write the array
*/
int *int2genotype( struct landscape h, int x , int *genotype);
int genotype2int( struct landscape h, int *genotype);

int int2nbmut( struct landscape h, int *x );

int isvalueinarray(int val, int *arr, int size);

/*
	init, free and print a landscape
	int *alleles is an array with the number of alleles per locus
	(not necessarily identical)
*/
void init_landscape( struct landscape *l, int nlocus, int *alleles);
void free_landscape( struct landscape *l );
void print_landscape( struct landscape *l );     /* give info about the landscape */
void output_landscape( struct landscape *l );     /* output only the landscape --can be read by ReadFile-- */


void exp_landscape( struct landscape *fl );
void log_landscape( struct landscape *fl );   /* change to exp/log scale all fitness values */
void remove_negative( struct landscape *fl, char opt_cut );   /* makez sure that no values are negative opt_cut==1, set to zero negative, otherwise shift all values to positive */


void set_mask(struct landscape *land,int *m);


#endif

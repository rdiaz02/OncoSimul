/*
	header for count_genotypes.c functions
	this include ths list definition
	
	Basic focus is a staring genotype -- count neighbors, retrieve them, etc.

*/


#ifndef _GENOTYPES_H_
#define _GENOTYPES_H_


struct list {
	long n;
	long *genotypes;
};

/*
	create a new empty liost
*/
struct list new_empty_list();


/*
	if inc>0, it increases size list
	if <0, it shrinks the list
	if n+inc equals 0, it frees the list
*/
void resize_list( int inc, struct list *l );

/*
	free the input list
*/
void free_list( struct list *l  );


/*
	return 1, 0 , -1 that corresponds to ( f2 larger ),  (both eq), ( f1 larger )
	the operator are modulated by the Fitness_Increase factor (set to 1 for absolute inc/dec/neutral)
*/
int compare_fitness( float f1, float f2, float Fitness_Increase );


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
int RetrieveNeighbors( struct landscape *fl, int g, char *visited_genotypes, struct list *li, char *opt_Fit, float FitnessRatio );



/*
	g is the starting genotype
	if opt_cumul=1, cumulate the number of visisted genotypes, otherwise report their current number (after n steps)
	MaxSteps, stop the counting after MaxSteps. If set 0, then stop when the paths have all reached a peak.
	it returns the number of genotypes reached after s steps (from 1 to MaxSteps)
*/
struct list count_genotypes( int g, struct landscape *fl, char opt_cumul, char opt_LessFit, int MaxSteps, float FitnessRatio );




int CountFitterNeighbors( struct landscape *fl, int g, float IncreaseRatio, int opt_strict );
int CountFitterGenotypes( struct landscape *fl, int g, float IncreaseRatio, int opt_strict );


int CountDefinedNeighbor( struct landscape *fl, int g);

#endif

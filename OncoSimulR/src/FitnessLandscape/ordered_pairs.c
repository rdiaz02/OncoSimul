/*
	This counts the ditribution of all ordered pairs
*/


#include "landscape.h"
#include "genotypes.h"
#include <stdlib.h>
#include <stdio.h>


int *genotype_diff( int g, int g2, struct landscape *fl, int *ndiff )
{

	int *geno1=int2genotype( *fl, g , NULL),   /* the array form of both genotypes */
	    *geno2=int2genotype( *fl, g2 , NULL);
	
	int l;        /* a counter for locus */
	
	int *diff;    /* the array for the different locus */
	

	if( g == g2 ){
		*ndiff = 0;
		return NULL;
	}
	
	for(*ndiff=0, l=0; l<fl->nlocus; l++){
		if( geno1[l] != geno2[l] )
			(*ndiff)++;
	}

	diff = malloc( (size_t) (*ndiff)*sizeof(int) );
	if(!diff)fprintf(stderr, "genotype_diff: cannot allocate diff, bye"), exit(3);

	for(*ndiff=0, l=0; l<fl->nlocus; l++)
		if( geno1[l] != geno2[l] ){
			diff[ (*ndiff)++ ] = l;
		}
	
	return diff;

}


/*
	fl is the fitness landscape
	g the starting genotype
	FitnessRatio the minimum increase between two fitnesses to be considered as significant
	
	Histo_SitePairs is a count of all pairs of successive sites that show a fitness increase
	
	This assumes that Histo_SitePairs is adequately set
*/
int Count_LocusPairs( struct landscape *fl, int g, float FitnessRatio, float *Histo_LocusPairs )
{


	int x, y;  /* a locus */

	struct list fitter_g,
	            fitter_g2;

	int l1,
	    l2;

	int *LocusDiff=NULL;
	int ndiff;


	int npairs=0;

	fitter_g = new_empty_list();
	fitter_g2 = new_empty_list();

	/*
	print_intgenotype( g, fl );
	printf("\n");
	*/
	RetrieveNeighbors( fl, g, NULL, &fitter_g, "f", FitnessRatio );

	
	for( x = 0; x < fitter_g.n  ; x++){
	
		/*
		printf(" ------ ");
		print_intgenotype(fitter_g.genotypes[x], fl);
		*/
	
		LocusDiff = genotype_diff( g, fitter_g.genotypes[x], fl, &ndiff );
		l1 = LocusDiff[0];
		free(LocusDiff);
		
		/*printf(" mut at locus %d (%d)\n", l1, ndiff);*/


		RetrieveNeighbors( fl, fitter_g.genotypes[x], NULL, &fitter_g2, "f", FitnessRatio );
		
		/* printf("%d better\n", fitter_g2.n );*/
		
		
		for( y=0; y < fitter_g2.n ; y++){
		
		
			/*printf(" ------ ------ ");
			print_intgenotype( fitter_g2.genotypes[y], fl );*/

			LocusDiff = genotype_diff( fitter_g.genotypes[x], fitter_g2.genotypes[y], fl, &ndiff );
			l2 = LocusDiff[0];
			free(LocusDiff);

			/*printf(" mut at locus %d\n", l2);*/

			Histo_LocusPairs[ fl->nlocus * l1 + l2  ] += 1.0 ; /* / (float)fitter_g2.n; */
			
			/*printf(" add %f to %d>%d\n", 1.0 / (float)fitter_g2.n, l1, l2);*/
			
		}
		
		if(fitter_g2.n)npairs++;
	
		free_list( &fitter_g2 );
	
	}
	
	
	free_list( &fitter_g );
	return npairs;

}

float *Count_AllLocusPairs( struct landscape *fl, float FitnessRatio  )
{

	
	int g;
	float *Histo_LocusPairs;
	int tot_pairs=0;
	
	
	Histo_LocusPairs = (float *)calloc( fl->nlocus*fl->nlocus, sizeof(float) );
	if(!Histo_LocusPairs)fprintf(stderr, "Count_AllLocusPairs: cannot allocate Histo_LocusPairs, bye\n"), exit(3);
	
	for( g=0; g< fl->ngenotypes; g++)
		tot_pairs += Count_LocusPairs( fl, g, FitnessRatio,  Histo_LocusPairs);
	

	/*for( p=0; p<fl->nlocus*fl->nlocus ; p++ )
		Histo_LocusPairs[ p ] /= (float)tot_pairs;*/

	return Histo_LocusPairs;

}

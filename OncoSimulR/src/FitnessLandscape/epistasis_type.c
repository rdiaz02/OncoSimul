#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "landscape.h"

#define EPSILON 1e-6

/*
	Return 0 for nonce, 1 for magnitude, 2 for sign et 3 for reciprocal sign
*/
short int get_sign_epistasis( int *geno00, int *geno11, struct landscape *fl )
{

	int g[4] = {-1, -1,-1, -1};              /* genotypes are stored in the bit order 00 01 10 11 */
	int i;
	char nmut=0;
	int tmp;
	
	float d_01_00,d_10_00,d_11_10,d_11_01 ;
	
	/*
	printf(" << ");
	print_genotype(geno00, fl->nlocus);
	printf(" vs ");
	print_genotype(geno11, fl->nlocus);
	printf(" >>\n");
	*/
	
	/*
		Extract the genotype indices of the two intermediates
	*/
	g[0] = genotype2int( *fl, geno00 );
	g[3] = genotype2int( *fl, geno11 );

	for( i=fl->nlocus-1; i>=0; i-- )
		if( geno00[i] != geno11[i] )
		{
			
			tmp = geno00[i];
			geno00[i] = geno11[i];                       /* put first mut of 11 in the background of 00 */
			g[1+nmut] = genotype2int( *fl, geno00 );     /* gets its index */
			geno00[i] = tmp;                             /* mutate back */
			
			nmut =1;
		}
	
	
	d_01_00 = fl->fitness[ g[ 1 ] ] - fl->fitness[ g[ 0 ] ];   /* 01 - 00 */
	d_10_00 = fl->fitness[ g[ 2 ] ] - fl->fitness[ g[ 0 ] ];   /* 10 - 00 */
	
	d_11_10 = fl->fitness[ g[ 3 ] ] - fl->fitness[ g[ 2 ] ];   /* 11 - 10 */
	d_11_01 = fl->fitness[ g[ 3 ] ] - fl->fitness[ g[ 1 ] ];   /* 11 - 01 */
	
	/*
		Check type (if any) epistasis --in logscale--
	*/
	
	if( fabs(d_01_00 - d_11_10 ) < EPSILON && fabs(d_10_00 - d_11_01 ) < EPSILON )
		return 0;
	
	if( d_01_00 * d_11_10 >=0 && d_10_00 * d_11_01 >=0)
		return 1;
	
	if( d_01_00 * d_11_10 <0 && d_10_00 * d_11_01 <0)
		return 3;
	
	return 2;
	
	
}




int *get_sign_epistasis_FL( struct landscape *fl ){


	int *epistasis;

	int g;
	int *geno=NULL, *geno2=NULL;
	
	int l1, l2;     /* locus to be mutated */
	int a1, a2;     /* alleles of each locus */
	
	// int flag=0;
	
	epistasis = calloc( (size_t)4 , (size_t) sizeof(int) );
	if(!epistasis)fprintf(stderr, "get_sign_epistasis_FL: cannot allocate epistasis, bye\n" ), exit(3);
	
	/*if(fl->log_scale == 0 ){
		log_landscape( fl );
		flag=1;
	}*/
	
	for(g=0;g<fl->ngenotypes; g++)
	{
	
		geno = int2genotype( *fl , g, geno );
		geno2 = int2genotype( *fl , g, geno2 );
		
		for(l1 = fl->nlocus-1; l1>=0; l1--){
		
			for( a1=geno[l1]+1; a1<fl->alleles[l1]; a1++ )
			{
				geno2[l1] = a1;
				
				for(l2=l1-1; l2>=0; l2--)
				{
					for( a2=geno[l2]+1; a2<fl->alleles[l2]; a2++ )
					{
						geno2[l2] = a2;
						epistasis[ get_sign_epistasis( geno, geno2, fl ) ]++;
					}
					geno2[l2] = geno[l2];
				}
			
			}
			
			geno2[l1] = geno[l1];
		}
	}

	/*
	if( flag==1 )
		exp_landscape( fl );
	*/

	return epistasis;

}

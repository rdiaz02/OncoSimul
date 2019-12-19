/*
	This file lists all functions to generate
	model FL
*/



#ifndef _MODELS_H_
#define _MODELS_H_



/* extern char verbose; */



struct model_opt {


	/*
		A fixed amount
	*/
	float fix;

	/*
		Multiplicative
	*/
	float mu_s;   /* mean of selective coef -- fitness is then 1+s -- */
	float sigma_s;
	float DimRet;  /*  coeficient that is multiplied for each mutation from the 0000 genotype */
	
	/*
		House of Cards
	*/
	float sigma_hoc;   /* variance for the noise */

	/*
		Kaufmann NK
	*/
	int kaufman_K;       /* number of interacting loci */
	int kaufman_rand;    /* chosen at random ? */

	
	
	/*
		Eggbox
	*/
	float mu_eggbox;   /* fitness are +/- e with stdev E */
	float sigma_eggbox;


	/*
		Ising
	*/
	float mu_ising;
	float sigma_ising;
	short circular_ising;
	
	/*
		Dininushing Return (if d<1) / Magnitude Epistasis (when d>1)
	*/
	
	/*
		Optimal product
	*/
	float mu_optimum;        /* the ideal production */
	float sigma_optimum;     /* how fitness decreases */
	float mu_prod;           /* the mean prod per locus */
	float sigma_prod;        /* the stdev of prod per locus */
};


struct model_opt init_model( void  );


void MixedModel( struct landscape *l, struct model_opt options, int opt_log );
void DiminishingReturn( struct landscape *l, float power, int opt_log );
void EggBox( struct landscape *l, float mu_e, float sigma_e  );
void Ising( struct landscape *l, float mu_c, float sigma_c, short opt_circ  );

/*
	Kaufman NK landscapes where N is the locus number
	and K is the number of interacting loci per loci (K \in [0, N-1])
	if random=0, the neighbors are chosen
	if random=1, chose randomly interacting loci (with no symetry).
*/
void Kaufman_NK( struct landscape *l,  int K, short random );

/*
	Each locus contribue multiplicatively, but the fitness is a gaussian to an optimum value
*/
void Optimum( struct landscape *l , double mean_p, double sigma_p, double mean_o, double sigma_o );

/*
	It generates a landscape where each genotype has a multiplicative contribution
	-- additive in logscale ; this corresponds to NK landscape with K=0 --

	s gives a fitness of 1+s in normal scale, and simply s in logscale
*/
void Multiplicative( struct landscape *l , double mu, double sigma, int  opt_single_fitness, float DimRet, char opt_log );

/*
	Generate a landscape with complete epistasis
	equivalent in regards of fitness ranks to Kaufman_NK with K=N-1
	--faster implementation though--
*/
void HouseOfCards( struct landscape *l, float sigma);





/*
	OLD MODELS
*/


/*
	Built a Ease spinglass
	
	I is the number of interacting loci
	c is the weight of incompatibility between 2 locus
	random if interacting loci are random
*/

void Spinglass(struct landscape *l,int I , float c, short random );




/*
	Incompatibility landscapes where N is the locus number
	and K is the number of interacting loci per loci (K \in [0, N-1])
	there is a cost with interacting loci have not the same allele number
	
	if random=0, the neighbors are chosen
	if random=1, chose randomly interacting loci (with no symetry).
	
	if ordered=0, the fitness cost is (1/2)^x, wfor(g=0;g<l->ngenotypes;g++){for(g=0;g<l->ngenotypes;g++){
	where x is the locus number
	
	if linear is 1, the first and the last locus cannot have interactions together
	
	if ordered=1
		--> incompatibilities are set to (1/2)^(i+1) with i the locus number

	if mean != 0
		--> all incompatibilities have a cost of mean if stdev == 0
		--> incompatibilities have a cost of Normal(mean, stdev) otherwise
	
	if mean == 0 && stdev == 0
		--> incompatibilities have a cost of Uniform[0,1]
	
*/
void Incompatibilities( struct landscape *l,  int K, short random, short ordered, short opt_linear, float mean, float stdev );



#endif

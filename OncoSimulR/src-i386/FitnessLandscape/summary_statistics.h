#ifndef _GAMMA_H_
#define _GAMMA_H_

#include "LinearAlgebra.h" // :: NEEDED HERE FOR STRUCT MATRIX DEFINITION ::

/*
	output stats
*/

void OutputSummaryStats( struct landscape * FL, float IncreaseRatio, char opt_log , char opt_short, char *infile, char *outfile );   /* this is a texte output */
 
char * outputstats( struct landscape * FL, float IncreaseRatio, char opt_log , char *FLname, char *CVSname);


/*
	Sinks and Peaks
*/
int numberPeaks(struct landscape *fl, float IncreaseRatio);
int numberSinks(struct landscape *fl, float IncreaseRatio);



/*
	CHAINS
	Function and structure used for chains

	In chains.c
*/

struct FromTo {

	long To;        /* only 1 way out */
	long *From;     /* but several origins From[0] contains the number
	                   of origins and From[1] to From[ From[0] ] the origins */
};

struct Chains {
	long nchains;
	long *depth;
	long *steps;
	long *origins;
};


int *compute_HistoBetterNeighbor( struct landscape * FL, float FitnessRatio );

void Connect_Fitness( struct FromTo * Ptr, struct landscape *fl, int g, float IncreaseRatio );

struct Chains Compute_Chains( struct landscape *fl ,  short opt_sum, float Increase );

long GetStartingPoints(  struct FromTo *Ptr, int g, struct landscape *fl );
long GetDepth(  struct FromTo *Ptr, int g,  struct landscape *fl );
long GetSteps(  struct FromTo *Ptr, int g,  struct landscape *fl );
void EraseChain( struct FromTo *Ptr, int g, struct landscape *fl  );

void free_chain(  struct Chains mychain );




/*
	Roughness / slope
	In summary statistics.c
*/

float compute_rs(  struct landscape *fl );

void compute_TransitionMatrix(struct landscape * fl, struct matrix * P, int * noOptima);
int compute_HammingDistance( int *geno1, int *geno2, int noLoci );

void compute_QMatrix(struct matrix *P, struct matrix *Q);
void compute_RMatrix(struct matrix *P, int noOptima, int * optimaIndex, struct matrix *R);
void compute_FundamentalMatrix(struct matrix *Q, struct matrix *N);
void compute_ExpectedNumberOfSteps(struct matrix *N, struct matrix *expectedNoSteps);
void compute_VarianceNumberOfSteps(struct matrix *N, struct matrix *expectedNoSteps, struct matrix *varianceNoSteps);
void compute_TransientProbs(struct matrix *N, struct matrix *transientProbs);
void compute_AbsorbingProbs(struct matrix *N, struct matrix *R, struct matrix *absorbingProbs);
int computeReachibility(struct matrix *transientProbs, struct matrix *B, int g);

/*
	In ordered_pairs.c
*/

int *genotype_diff( int g, int g2, struct landscape *fl, int *ndiff );
void Count_LocusPairs( struct landscape *fl, int g, float FitnessRatio, float *Histo_LocusPairs );
float *Count_AllLocusPairs( struct landscape *fl, float FitnessRatio  );



/*

	In gamma.c
*/

unsigned long myPow(unsigned long x, int p);

int gmod(struct landscape h, int g, int i);
int gmodMultiAllele(struct landscape h, int g, int locus, int allele);

double deltaf(struct landscape h, int g, int j, double tolerance);
double deltafMultiAllele(struct landscape h, int g, int locus, int allele, double tolerance);

//void GammaMultiAllele_ij(struct landscape fl, struct matrix * Nominator, struct matrix * Denominator, double tolerance);
void GammaMultiAllele_AiBiAjBj(struct landscape fl, struct matrix * NominatorAiBiAjBj, struct matrix * DenominatorAiBiAjBj, struct matrix * NominatorAiBi, struct matrix * DenominatorAiBi, struct matrix * NominatorAjBj, struct matrix * DenominatorAjBj, struct matrix * Nominatorij, struct matrix * Denominatorij, double tolerance);
void GammaMultiAllele_Ai_Bj(struct landscape fl, int noAlleles, struct matrix * NominatorAiBi, struct matrix * DenominatorAiBi, struct matrix * NominatorAjBj, struct matrix * DenominatorAjBj, struct matrix * NominatorAi, struct matrix * DenominatorAi, struct matrix * NominatorBj, struct matrix * DenominatorBj, double tolerance);
void GammaMultiAllele_AiBj(struct landscape fl, int noAlleles, struct matrix * NominatorAiBiAjBj, struct matrix * DenominatorAiBiAjBj, struct matrix * NominatorAiBj, struct matrix * DenominatorAiBj, double tolerance);
void GammaMultiAllele_ij(struct landscape fl, int noAlleles, struct matrix * NominatorAi, struct matrix * DenominatorAi, struct matrix * NominatorBj, struct matrix * DenominatorBj, struct matrix * Nominatori, struct matrix * Denominatori, struct matrix * Nominatorj, struct matrix * Denominatorj, double tolerance);

double GammaDistance(struct landscape fl, int distance, int tolerance);

double gamma_den(struct landscape h,double tolerance);
double gamma_den_i(struct landscape h, int i,double tolerance);
double gamma_den_j(struct landscape h, int j,double tolerance);

double gamma_num(struct landscape h, int g,double tolerance);
double gamma_num_i(struct landscape h, int g, int i,double tolerance);
double gamma_num_j(struct landscape h, int g, int j,double tolerance);
double gamma_num_ij(struct landscape h, int g, int i, int j,double tolerance);

double gamma_num_global(struct landscape h,double tolerance);
double gamma_num_i_global(struct landscape h, int i,double tolerance);
double gamma_num_j_global(struct landscape h, int j,double tolerance);
double gamma_num_ij_global(struct landscape h, int i, int j,double tolerance);
//distance greater than 1:

double gamma_num_dist_j_global(struct landscape h, int d, int j,double tolerance);
double gamma_num_dist_global(struct landscape h, int d,double tolerance);

double gamma_global(struct landscape h,double tolerance);
double gamma_i_global(struct landscape h, int i,double tolerance);
double gamma_j_global(struct landscape h, int j,double tolerance);
double gamma_ij_global(struct landscape h, int i, int j,double tolerance);
//distance greater than 1:

double gamma_dist_j_global(struct landscape h, int d, int j,double tolerance);
double gamma_dist_global(struct landscape h, int d,double tolerance);



/*
	In epistasis_type.c
*/
/*
	Return 0 for magnitude, 1 for sign et 2 for reciprocal sign
*/
short int get_sign_epistasis( int *geno00, int *geno11, struct landscape *fl );
int *get_sign_epistasis_FL( struct landscape *fl );

/*
	In decomposition.c
*/

void FourierTransform(struct landscape *FL, struct matrix * fourierCoefficients, double * fourierCoefficient_total, double * fourierCoefficient_epi, int * max_f);
int computeSubsetSize(struct landscape *FL, int locus, int currentSize, int maxSize);
void createSubsets(struct landscape *FL, int locus, int allele, int setIndex, int currentSize, int maxSize, int **alleleArray, int *** subsets);

#endif




/*

	Where we store all summary statistics to describe the F.L.

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <strings.h>
#include <string.h>

#include "landscape.h"
#include "genotypes.h"
#include "summary_statistics.h"
#include "calculus.h"

#define EPSILON 1e-10
#define THRESHOLD 1e-6

/*
	Return number of peaks in the given landscape
*/
int numberPeaks(struct landscape * fl, float IncreaseRatio)
{
	int g=0,     /* current genoytpe */
	    nbp=0;   /* number of peaks */
	
	
	for (g = 0; g < fl->ngenotypes; g++)
		if( CountFitterNeighbors(fl, g, IncreaseRatio, 1) == 0 )
			nbp++;

	return nbp;
}


/*
	Return number of sinks in the given landscape
*/
int numberSinks(struct landscape * fl, float IncreaseRatio)
{
	int             g = 0, nbs = 0;
	for (g = 0; g < fl->ngenotypes; g++)
		if (CountFitterNeighbors(fl, g, IncreaseRatio, 1) == fl->neighbors)
			nbs++;
	return nbs;
}


void compute_TransitionMatrix(struct landscape * fl, struct matrix *P, int * noOptima)
{
	
	int *geno1 = NULL;
	int *geno2 = NULL;
	int HD;
	float normalizeConst;
	float selco;
	int g1, g2;
	
	
	for(g1 = 0; g1 < P->r; g1++)
	{
		geno1 = int2genotype(*fl, g1, geno1);
		normalizeConst = 0;

		for(g2 = 0; g2 < P->r; g2++)
		{
			geno2 = int2genotype(*fl, g2, geno2);
			HD = compute_HammingDistance(geno1, geno2, fl->nlocus);
			if (HD == 1)
			{
				// :: WHICH SELCO DEF TO USE HERE :: DEPENDS ON SCALE :: //
				selco = fl->fitness[g2] - fl->fitness[g1];

				/*
				selco = fl->fitness[g2] / fl->fitness[g1] - 1.;
				*/
				
				if (selco > 0)
				{
					P->val[g1][g2] = 1. - exp(-2 * selco);
					normalizeConst += 1. - exp(-2 * selco);
				}
				else
				{
					P->val[g1][g2] = 0.;
				}
			}
			else
			{
				P->val[g1][g2] = 0.;
			}
		}
		if (normalizeConst == 0)
		{
			P->val[g1][g1] = 1.;
			(*noOptima)++;
			
		}
		else
		{
			for(g2 = 0; g2 < P->c; g2++)
			{
				P->val[g1][g2] /= normalizeConst;
			}
		}
	}
	
	free(geno1);
	free(geno2);
}

// :: THIS IS BASICALLY A SUBPART OF 'int *genotype_diff( int g, int g2, struct landscape *fl, int *ndiff )' :: //
// :: COULD POTENTIALLY BE INTEGRATED/DROPPED IF NEED BE :: //

int compute_HammingDistance(int *geno1, int *geno2, int noLoci)
{
	int HD = 0, i;

	for(i = 0; i < noLoci; i++)
	{
		if (geno1[i] != geno2[i])
		{
			HD++;
		}
	}

	return(HD);
}


void compute_QMatrix(struct matrix *P, struct matrix *Q)
{

	int i, j, index1, index2;

	index1 = 0;
	for(i = 0; i < P->r; i++)
	{
		index2 = 0;
		for(j = 0; j < P->r; j++)
		{
			if ( (P->val[i][i] != 1.) && (P->val[j][j] != 1.) )
			{
				Q->val[index1][index2] = P->val[i][j];
				index2++;
			}
		}
		if (P->val[i][i] != 1.)
		{
			index1++;
		}
	}
}

void compute_RMatrix(struct matrix *P, int noOptima, int * optimaIndex, struct matrix *R)
{
	int i, j, index1, index2;
	
	index1 = 0;
	for(i = 0; (i < P->r) || ( index1 == (P->r - noOptima - 1)); i++)
	{
		if (P->val[i][i] != 1.)
		{
			index2 = 0;
			for(j = 0; (j < P->c) || (index2 == (noOptima - 1)); j++)
			{
				
				//printf("\t%d", j);
				if ( (P->val[i][i] != 1.) && (P->val[j][j] == 1.) )
				{
					R->val[index1][index2] = P->val[i][j];
					optimaIndex[index2] = j;
					index2++;
				}
			}
			index1++;
		}
	}
	
}


void compute_FundamentalMatrix(struct matrix *Q, struct matrix *N)
{
	struct matrix I, temp;
	
	I = MemMat(Q->r, Q->c);
	IdentMat(&I);
	temp = MatrixDiff(&I, Q, 0);
	inverseMatrix(&temp, N);

	FreeMat(&I);
	FreeMat(&temp);
}

void compute_ExpectedNumberOfSteps(struct matrix *N, struct matrix *expectedNoSteps)
{
	struct matrix oneVector;
	int i;

	oneVector = MemMat(N->r, 1);
	for(i = 0; i < N->r; i++)
	{
		oneVector.val[i][0] = 1.;
	}

	MatrixProductRef(N, &oneVector, expectedNoSteps, 0);
	FreeMat(&oneVector);
}

void compute_VarianceNumberOfSteps(struct matrix *N, struct matrix *expectedNoSteps, struct matrix *varianceNoSteps)
{
	struct matrix I, expectedNoStepsSquared, temp1, temp2, temp3;
	int i;
	
	I = MemMat(N->r, N->c);
	IdentMat(&I);
	expectedNoStepsSquared = MatrixHadarmardProduct(expectedNoSteps, expectedNoSteps, 0);
	temp1 = MemMat(N->r, N->c);
	scaleMatrixRef(N, &temp1, 2., 0);
	temp2 = MemMat(N->r, N->c);
	MatrixDiffRef(&temp1, &I, &temp2, 0);
	temp3 = MatrixProduct(&temp2, expectedNoSteps, 0);
	
	MatrixDiffRef(&temp3, &expectedNoStepsSquared, varianceNoSteps, 0);

	for(i = 0; i < N->c; i++)
	{
		if(varianceNoSteps->val[i][0] < 0)
		{
			varianceNoSteps->val[i][0] = 0.;
		}
	}
	
	FreeMat(&expectedNoStepsSquared);
	FreeMat(&temp2);
	FreeMat(&temp3);
}

void compute_TransientProbs(struct matrix *N, struct matrix *transientProbs)
{
	struct matrix temp1, temp2, temp3, I;
	I = MemMat(N->r, N->c);
	IdentMat(&I);
	
	temp1 = MatrixDiff(N, &I, 0);
	temp2 = CreateMatrixDiagMatrix(N, 0);
	temp3 = Invert_MatrixDiag(&temp2, 0);
	
	MatrixProductRef(&temp1, &temp3, transientProbs, 1);

	FreeMat(&temp2);
	FreeMat(&I);
}

int computeReachibility(struct matrix *transientProbs, struct matrix *B, int g)
{
	int count = 0, i;
	for(i = 0; i < transientProbs->c; i++)
	{
		if(transientProbs->val[g][i] > THRESHOLD)
		{
			count++;
		}
	}
	
	for(i = 0; i < B->c; i++)
	{
		if(B->val[g][i] > THRESHOLD)
		{
			count++;
		}
	}
	
	return (count);
}


void compute_AbsorbingProbs(struct matrix *N, struct matrix *R, struct matrix *absorbingProbs)
{
	MatrixProductRef(N, R, absorbingProbs, 0);
}



float compute_rs(struct landscape * fl )
{

	/*
		Calculate roughness-to-slope ratio as described in Aita et al. 2001 Prot. Engineer.
	*/
	
	struct matrix   X, Y;
	struct matrix   beta;

	float           roughness;
	float           *slope = (float *)malloc( (size_t) fl->nlocus*sizeof(float) );
	float			meanSlope = 0;
	
	int             g;
	int            *geno = NULL;

	float           ftheo = 0;

	int             i,j, cAlleles;

	/*
		First build the X matrix and the Y vector
	*/

	X = MemMat(fl->ngenotypes, fl->nalleles - fl->nlocus + 1);	/* this one has all genotypes */
	
	Y = MemMat(fl->ngenotypes, 1);	/* this one will contain the fitness */

	beta = MemMat(fl->nlocus + 1, 1);

	for (g = 0; g < fl->ngenotypes; g++)
	{
		geno = int2genotype(*fl, g, geno);
		cAlleles = 0;
		for (j = 0; j < fl->nlocus; j++)
		{
			for (i = 1; i < fl->alleles[j]; i++)
			{
				if(geno[j] == i)
				{
					X.val[g][cAlleles] = 1.;
				}
				else
				{
					X.val[g][cAlleles] = 0.;
				}
				
				cAlleles++;
			}
		}
		X.val[g][cAlleles] = 1.0;	/* the basal fitness */
		Y.val[g][0] = (float) fl->fitness[g];

	}

	/*
		in beta[i], the predicted increase of locus i
		in beta[nlocus], is stored rthe basal fitness
	*/
	
	
	beta = LeastSquare(&X, &Y);
	
	/*
		Compute slope
	*/
	cAlleles = 0;
	for (j = 0; j < fl->nlocus; j++)
	{
		slope[j] = 0.0;
		for (i = 1; i < fl->alleles[j]; i++)
		{
			slope[j] += (float) fabs(beta.val[cAlleles][0]);
			cAlleles++;
		}
		slope[j] /= ((float) fl->alleles[j] - 1.);
		meanSlope += (float) fabs(slope[j]);
	}

	meanSlope /= (float) fl->nlocus;
	
	/*
		Compute noise, aka roughness, aka residuals
	*/
	roughness = 0;

	for (g = 0; g < fl->ngenotypes; g++)
	{

		geno = int2genotype(*fl, g, geno);
		cAlleles = 0;
		ftheo = beta.val[fl->nalleles - fl->nlocus][0];
		
		for (j = 0; j < fl->nlocus; j++)
		{
			for (i = 1; i < fl->alleles[j]; i++)
			{
				if(geno[j] == i)
				{
					ftheo += beta.val[cAlleles][0];
				}
				cAlleles++;
			}
		}

		roughness += pow(Y.val[g][0] - ftheo, 2.0);

	}
	

	
	roughness /= fl->ngenotypes;
	roughness = sqrt(roughness);
	
	free(geno);
	FreeMat(&X);
	FreeMat(&Y);
	FreeMat(&beta);
	free(slope);

	return (roughness / meanSlope);
}







/*
	
	outputs some nice values about landscape
	
	deprecated for command line version
	
*/

/*char * outputstats(struct landscape * land, int web)
{
	int             i;
	int             p;
	int             t = 1;
	char            endofline[10];
	char           *lres = NULL;
	double         *flux;
	int            *m;
	int             n = land->ngenotypes;

	flux = malloc(sizeof(double) * n * n);
	m = malloc(sizeof(int) * n * n);
	
	if (web == 1)
		strcpy(endofline, "<BR>\n");
	else
		strcpy(endofline, "\n");

	lres = malloc(sizeof(char) * 10024);
	
	
	
	if (web == 0)
		printf("Some statistics%s", endofline);
	else
		sprintf(lres, "Some statistics%s", endofline);
		
	p = numberPeaks(land, 1);
	
	if (web == 0)
		printf("Nbr of peaks=%d%s", p, endofline);
	else
		sprintf(lres, "%sNbr of peaks=%d%s", lres, p, endofline);
		
	if (web == 0)
		printf("Nbr of sinks=%d%s", numberSinks(land, 1), endofline);
	else
		sprintf(lres, "%sNbr of sinks=%d%s", lres, numberSinks(land, 1), endofline);
		
	if (web == 0)
		printf("Summary of gamma (only for 2 allels):%s", endofline);
	else
		sprintf(lres, "%sSummary of gamma (only for 2 allels):%s", lres, endofline);
		
	
	for (i = 0; i < land->nlocus; i++)
		if (land->alleles[i] != 2) {
			t = 0;
			break;
		}
		
		
	if (t == 1) {

		for (i = 1; i < land->nlocus; i++) {

			if (web == 0)
				printf(" gamma global (dist %d) = %f %s", i, gamma_dist_global(*land, i, -1), endofline);
			else
				sprintf(lres, "%s gamma global (dist %d) = %f %s", lres, i, gamma_dist_global(*land, i, -1), endofline);
		}

	}
	

	
	free(flux);
	free(m);
	free(lres);
	
	if (web == 1)
		return (lres);
	else
		return (NULL);
}
*/




char * outputstats( struct landscape * FL, float IncreaseRatio, char opt_log , char *FLname,char *CSVfile)
{
	int i, j, g;
	
	int npeaks;                       /* number of peaks and sinks */
	int nsinks;
	int noOptima = 0;
	double meanExpectedNumberOfStepsFL = 0;
	double meanVarianceNumberOfStepsFL = 0;
	double meanReachabilityFL = 0;
	double meanFitterGenotypesFL = 0;
	
	FILE *fcsv; 
	
	struct matrix P = MemMat(FL->ngenotypes, FL->ngenotypes);
	compute_TransitionMatrix(FL, &P, &noOptima);					/* Transition matrix for entire fitness landscape*/
	
	double meanReachOptFL[noOptima];
	for(i = 0; i < noOptima; i++)
	{
		meanReachOptFL[i] = 0;
	}

	struct matrix Q = MemMat(P.r - noOptima, P.c - noOptima);
	compute_QMatrix(&P, &Q);										/* Transition matrix between transitive states (non-optima)*/
	
	int * optimaIndex = (int *)malloc( (size_t) noOptima*sizeof(int) );
	
	struct matrix R = MemMat(P.r - noOptima, noOptima);
	compute_RMatrix(&P, noOptima, optimaIndex, &R);					/* Transition matrix between transitive states -> optima; also returns index of optima */
	
	struct matrix N	= MemMat(P.r - noOptima, P.c - noOptima);
	compute_FundamentalMatrix(&Q, &N);								/* Fundamental matrix N = (I-Q)^-1 */
	
	struct matrix expectedNoSteps = MemMat(N.r, 1);
	compute_ExpectedNumberOfSteps(&N, &expectedNoSteps);			/* Expected number of steps (until being absorbed/ any optimum is reached) t = N 1 */
	
	struct matrix varianceNoSteps = MemMat(N.r, 1);
	compute_VarianceNumberOfSteps(&N, &expectedNoSteps, &varianceNoSteps);	/*Variance in number of steps (until being absorbed/ any optimum is reached) = (2N-I)t - t'*'t ('*': Hadarmard product) */
	
	struct matrix transientProbs = MemMat(N.r, N.c);
	compute_TransientProbs(&N, &transientProbs);					/* Transient probabilities to reach any non-absorbing state from any non-absorbing state (eventually) */
	
	struct matrix B = MemMat(P.r - noOptima, noOptima);
	compute_AbsorbingProbs(&N, &R, &B);								/* Absorbing probabilities to reach any optimum from any non-absorbing state (eventually)*/
	
	struct matrix reachibility = MemMat(Q.r, 1);
	struct matrix fitterGenotypes = MemMat(Q.r, 1);
	
	i = 0;
	for (g = 0; g < FL->ngenotypes; g++)
	{
		if (P.val[g][g] != 1.)
		{
			reachibility.val[i][0] = computeReachibility(&transientProbs, &B, i);	/* Check if given transient state can be reached from a given starting genotype, i.e., when their transient probability is > 0 */
			fitterGenotypes.val[i][0] = CountFitterGenotypes(FL, g, 1., 1);			/* Count the number of mutants with higher fitness than that of a given starting genotype */
			
			meanExpectedNumberOfStepsFL += expectedNoSteps.val[i][0];
			meanVarianceNumberOfStepsFL += varianceNoSteps.val[i][0];
			meanReachabilityFL += reachibility.val[i][0];
			meanFitterGenotypesFL += fitterGenotypes.val[i][0];
			
			for (j = 0; j < noOptima; j++)
			{
				meanReachOptFL[j] += B.val[i][j];
			}
			i++;
		}
	}
	
	meanExpectedNumberOfStepsFL /= (0.0 + FL->ngenotypes - noOptima);
	meanVarianceNumberOfStepsFL /= (0.0 + FL->ngenotypes - noOptima);
	meanReachabilityFL /= (0.0 + FL->ngenotypes - noOptima);
	meanFitterGenotypesFL /= (0.0 + FL->ngenotypes - noOptima);
	
	for (j = 0; j < noOptima; j++)
	{
		meanReachOptFL[j] /= (0.0 + FL->ngenotypes - noOptima);
	}
	
	int sizeSubSets = 0;
	int noAlleles = 0;
	for (j = 0; j < FL->nlocus; j++)
	{
		sizeSubSets += N_choose_n(FL->alleles[j], 2);
		noAlleles += FL->alleles[j];
	}
	
	
	struct matrix GammaMultiAllele_AiBiAjBj_Nominator = MemMat(sizeSubSets, sizeSubSets);
	struct matrix GammaMultiAllele_AiBiAjBj_Denominator = MemMat(sizeSubSets, sizeSubSets);
	struct matrix GammaMultiAllele_AiBiAjBj_Matrix = MemMat(sizeSubSets, sizeSubSets);
	
	struct matrix GammaMultiAllele_AiBi_Nominator = MemMat(sizeSubSets, 1);
	struct matrix GammaMultiAllele_AiBi_Denominator = MemMat(sizeSubSets, 1);
	struct matrix GammaMultiAllele_AiBi_Matrix = MemMat(sizeSubSets, 1);
	
	struct matrix GammaMultiAllele_AjBj_Nominator = MemMat(sizeSubSets, 1);
	struct matrix GammaMultiAllele_AjBj_Denominator = MemMat(sizeSubSets, 1);
	struct matrix GammaMultiAllele_AjBj_Matrix = MemMat(sizeSubSets, 1);
	
	struct matrix GammaMultiAllele_AiBj_Nominator = MemMat(noAlleles, noAlleles);
	struct matrix GammaMultiAllele_AiBj_Denominator = MemMat(noAlleles, noAlleles);
	struct matrix GammaMultiAllele_AiBj_Matrix = MemMat(noAlleles, noAlleles);
	
	struct matrix GammaMultiAllele_ij_Nominator = MemMat(FL->nlocus, FL->nlocus);
	struct matrix GammaMultiAllele_ij_Denominator = MemMat(FL->nlocus, FL->nlocus);
	struct matrix GammaMultiAllele_ij_Matrix = MemMat(FL->nlocus, FL->nlocus);
	
	GammaMultiAllele_AiBiAjBj( *FL, &GammaMultiAllele_AiBiAjBj_Nominator, &GammaMultiAllele_AiBiAjBj_Denominator, &GammaMultiAllele_AiBi_Nominator, &GammaMultiAllele_AiBi_Denominator, &GammaMultiAllele_AjBj_Nominator, &GammaMultiAllele_AjBj_Denominator, &GammaMultiAllele_ij_Nominator, &GammaMultiAllele_ij_Denominator,-1);
	for (i = 0; i < sizeSubSets; i++)
	{
		for (j = 0; j < sizeSubSets; j++)
		{
			GammaMultiAllele_AiBiAjBj_Matrix.val[i][j] = GammaMultiAllele_AiBiAjBj_Nominator.val[i][j] / GammaMultiAllele_AiBiAjBj_Denominator.val[i][j];
		}
		GammaMultiAllele_AiBi_Matrix.val[i][0] = GammaMultiAllele_AiBi_Nominator.val[i][0] / GammaMultiAllele_AiBi_Denominator.val[i][0];
		GammaMultiAllele_AjBj_Matrix.val[i][0] = GammaMultiAllele_AjBj_Nominator.val[i][0] / GammaMultiAllele_AjBj_Denominator.val[i][0];
	}
	
	struct matrix GammaMultiAllele_Ai_Nominator = MemMat(noAlleles, 1);
	struct matrix GammaMultiAllele_Ai_Denominator = MemMat(noAlleles, 1);
	struct matrix GammaMultiAllele_Ai_Matrix = MemMat(noAlleles, 1);
	
	struct matrix GammaMultiAllele_Bj_Nominator = MemMat(noAlleles, 1);
	struct matrix GammaMultiAllele_Bj_Denominator = MemMat(noAlleles, 1);
	struct matrix GammaMultiAllele_Bj_Matrix = MemMat(noAlleles, 1);
	
	GammaMultiAllele_Ai_Bj( *FL, noAlleles, &GammaMultiAllele_AiBi_Nominator, &GammaMultiAllele_AiBi_Denominator, &GammaMultiAllele_AjBj_Nominator, &GammaMultiAllele_AjBj_Denominator, &GammaMultiAllele_Ai_Nominator, &GammaMultiAllele_Ai_Denominator, &GammaMultiAllele_Bj_Nominator, &GammaMultiAllele_Bj_Denominator, -1);
	GammaMultiAllele_AiBj( *FL, noAlleles, &GammaMultiAllele_AiBiAjBj_Nominator, &GammaMultiAllele_AiBiAjBj_Denominator, &GammaMultiAllele_AiBj_Nominator, &GammaMultiAllele_AiBj_Denominator, -1);
	for (i = 0; i < noAlleles; i++)
	{
		for (j = 0; j < noAlleles; j++)
		{
			GammaMultiAllele_AiBj_Matrix.val[i][j] = GammaMultiAllele_AiBj_Nominator.val[i][j] / GammaMultiAllele_AiBj_Denominator.val[i][j];
		}
		GammaMultiAllele_Ai_Matrix.val[i][0] = GammaMultiAllele_Ai_Nominator.val[i][0] / GammaMultiAllele_Ai_Denominator.val[i][0];
		GammaMultiAllele_Bj_Matrix.val[i][0] = GammaMultiAllele_Bj_Nominator.val[i][0] / GammaMultiAllele_Bj_Denominator.val[i][0];
	}
	
	struct matrix GammaMultiAllele_i_Nominator = MemMat(FL->nlocus, 1);
	struct matrix GammaMultiAllele_i_Denominator = MemMat(FL->nlocus, 1);
	struct matrix GammaMultiAllele_i_Matrix = MemMat(FL->nlocus, 1);
	
	struct matrix GammaMultiAllele_j_Nominator = MemMat(FL->nlocus, 1);
	struct matrix GammaMultiAllele_j_Denominator = MemMat(FL->nlocus, 1);
	struct matrix GammaMultiAllele_j_Matrix = MemMat(FL->nlocus, 1);
	
	
	GammaMultiAllele_ij( *FL, noAlleles, &GammaMultiAllele_Ai_Nominator, &GammaMultiAllele_Ai_Denominator, &GammaMultiAllele_Bj_Nominator, &GammaMultiAllele_Bj_Denominator, &GammaMultiAllele_i_Nominator, &GammaMultiAllele_i_Denominator, &GammaMultiAllele_j_Nominator, &GammaMultiAllele_j_Denominator, -1);
	for (i = 0; i < FL->nlocus; i++)
	{
		for (j = 0; j < FL->nlocus; j++)
		{
			GammaMultiAllele_ij_Matrix.val[i][j] = GammaMultiAllele_ij_Nominator.val[i][j] / GammaMultiAllele_ij_Denominator.val[i][j];
		}
		GammaMultiAllele_i_Matrix.val[i][0] = GammaMultiAllele_i_Nominator.val[i][0] / GammaMultiAllele_i_Denominator.val[i][0];
		GammaMultiAllele_j_Matrix.val[i][0] = GammaMultiAllele_j_Nominator.val[i][0] / GammaMultiAllele_j_Denominator.val[i][0];
	}
	
	//float gamma, gamma_star;                       /* general measure of FL curvature */

	double * gamma = (double *)malloc( (size_t) FL->nlocus*sizeof(double) );
	double * gamma_star = (double *)malloc( (size_t) FL->nlocus*sizeof(double) );
	
	int *HistoOutDegree=NULL;
	int nbBl=1;

	float rs;
	
	int nbBloks=10024;

	char *f;

	float none;
	int *e, total_e;
	
	
	struct Chains mychain;            /* structure to store chain results */
	long nsteps=0,
		depthmax =0,

	    norigins=0;
	
	f=malloc(sizeof(char )*(nbBloks));
	nbBl++;
	
	if(opt_log)
	{
		log_landscape( FL );
	}
		
	npeaks = numberPeaks( FL, IncreaseRatio );
	nsinks = numberSinks( FL, IncreaseRatio );
	
	gamma[0] = 1;
	gamma_star[0] = 1;
	for (i = 1; i < FL->nlocus; i++)
	{
		gamma[i] = GammaDistance(*FL, i, -1);
		gamma_star[i] = GammaDistance(*FL, i, 0);
	}

	
	HistoOutDegree = compute_HistoBetterNeighbor( FL, IncreaseRatio );
	
	mychain = Compute_Chains( FL, 0, IncreaseRatio  );
	
	rs = compute_rs(  FL );

	e = get_sign_epistasis_FL( FL );
	total_e=e[0]+e[1]+e[2]+e[3];

	

	// :: FOURIER DECOMPOSITION ::
	
	double fourierCoefficient_total = 0;
	double fourierCoefficient_epi = 0;
	
	struct matrix fourierCoefficients = MemMat(FL->nlocus + 1, 1);
	int max_f = 1;
	
	for (i = 0; i < FL->nlocus + 1; i++)
	{
		fourierCoefficients.val[i][0] = 0;
	}
	FourierTransform(FL, &fourierCoefficients, &fourierCoefficient_total, &fourierCoefficient_epi, &max_f);

	for( i=0 ; i < mychain.nchains ; i++ ){
		nsteps += mychain.steps[i];
		norigins += mychain.origins[i];
		depthmax = (mychain.depth[i]<depthmax)?depthmax:mychain.depth[i];

	}

	/****************************
		::  CSV output   ::
	 ***************************
	*/
	 	fcsv=fopen(CSVfile,"w");
	 	if (fcsv==NULL) printf("Pb with file %s<BR>",CSVfile),exit(1);
	 
	 	fprintf(fcsv,"genotypes, peaks, sinks, r/s, gamma, gamma*, none, magnitude, sign, recip. sign, chains, step, origins, max_depth, mean no. Steps, mean reachibility, mean no. fitter genotypes, nboptimum");
	 	for (i=1; i<=noOptima; i++)
			{
				fprintf(fcsv,", optimum_%d_index, optimum_%d_mean, optimum_%d_absorb, optimum_%d_proba",i,i,i,i);
			}
		fprintf(fcsv,",");		
	 	if (fourierCoefficient_total!=0)
			{
				
				for(i=1;i<FL->nlocus+1; i++)
					fprintf(fcsv,", fraction_%d",i);
	 		}
	 	fprintf(fcsv,"\n");// last item
	 	
	 	fprintf(fcsv,"%d,%d,%d,",FL->ngenotypes,npeaks, nsinks);
	 	fprintf(fcsv,"%.3f, %.3f, %.3f,",rs, gamma[1], gamma_star[1]);
	 	
	 	none=1-( (e[1]/(0.0+total_e))+(e[2]/(0.0+total_e))+(e[3]/(0.0+total_e)));
	 	
	 	fprintf(fcsv,"%.3f, %.3f, %.3f, %.3f,",none,e[1]/(0.0+total_e),e[2]/(0.0+total_e),e[3]/(0.0+total_e));
		fprintf(fcsv,"%ld, %ld, %ld, %ld,",mychain.nchains,nsteps,norigins,depthmax);
		fprintf(fcsv,"%.2f, %.2f, %.2f,",meanExpectedNumberOfStepsFL,meanReachabilityFL,meanFitterGenotypesFL);
		
		fprintf(fcsv,"%d,",noOptima);
		for (i=0; i<noOptima; i++)
			fprintf(fcsv,"%d, , , %.2f,",optimaIndex[i],meanReachOptFL[i]);//mean and abs are missing
		
		if (fourierCoefficient_total!=0)
			{
				
				for(i=1;i<FL->nlocus+1; i++)
					fprintf(fcsv,",%f",fourierCoefficients.val[i][0]/fourierCoefficient_total);
	 		}	
		
		
	 	fprintf(fcsv,"\n");
	 	fclose(fcsv);
	/****************************
		::   Brief Output   ::
	 ****************************/
		//https://stackoverflow.com/a/2674354
		// https://github.com/dyne/frei0r/issues/22
		// https://stackoverflow.com/questions/59708584/how-to-sprintf-to-write-a-string-without-warnings-about-restrictions
		sprintf(f,"&nbsp;&nbsp;General:<BR>&nbsp;&nbsp;&nbsp;#gentotype:%d&nbsp;&nbsp;&nbsp;<BR>&nbsp;&nbsp;&nbsp;#peaks:%d<BR>&nbsp;&nbsp;&nbsp;#sinks:%d<BR><BR>",FL->ngenotypes,npeaks, nsinks);
		// In what follows, I turn all things like the following, to the one below. I.e., + strlen(f), delete %s, and delete f after the ""
		// sprintf(f,"%s&nbsp;&nbsp;Amount of Epistasis:&nbsp;&nbsp;<BR>&nbsp;&nbsp;&nbsp;r/s:%.3f<BR>&nbsp;&nbsp;&nbsp;gamma:%.3f<BR>&nbsp;&nbsp;&nbsp;gamma*:%.3f<BR><BR>",f,rs, gamma[1], gamma_star[1]);
		sprintf(f + strlen(f),"&nbsp;&nbsp;Amount of Epistasis:&nbsp;&nbsp;<BR>&nbsp;&nbsp;&nbsp;r/s:%.3f<BR>&nbsp;&nbsp;&nbsp;gamma:%.3f<BR>&nbsp;&nbsp;&nbsp;gamma*:%.3f<BR><BR>",rs, gamma[1], gamma_star[1]);
		none=1-( (e[1]/(0.0+total_e))+(e[2]/(0.0+total_e))+(e[3]/(0.0+total_e)));
		sprintf(f + strlen(f),"&nbsp;Type of Epistasis:<BR>&nbsp;&nbsp;&nbsp;none:%.3f<BR>&nbsp;&nbsp;&nbsp;magnitude:%.3f<BR>&nbsp;&nbsp;&nbsp;sign:%.3f<BR>&nbsp;&nbsp;&nbsp;recipr. sign:%.3f<BR><BR>",none,e[1]/(0.0+total_e),e[2]/(0.0+total_e),e[3]/(0.0+total_e));
		sprintf(f + strlen(f),"&nbsp;Chains:<BR>&nbsp;&nbsp;&nbsp;chains:%ld<BR>&nbsp;&nbsp;&nbsp;steps:%ld<BR>&nbsp;&nbsp;&nbsp;origins:%ld<BR>&nbsp;&nbsp;&nbsp;max_depth:%ld<BR><BR>",mychain.nchains,nsteps,norigins,depthmax);
		
		// :: ADAPTIVE WALK BRIEF OUTPUT ::
		sprintf(f + strlen(f),"&nbsp;Adaptive Walks:<BR>&nbsp;&nbsp;&nbsp;mean no. Steps:%.2f<BR>&nbsp;&nbsp;&nbsp;mean reachibility:%.2f<BR>&nbsp;&nbsp;&nbsp;mean no. fitter genotypes:%.2f<BR><BR>",meanExpectedNumberOfStepsFL,meanReachabilityFL,meanFitterGenotypesFL);
		sprintf(f + strlen(f),"&nbsp;Optima:<BR>");
		for (i=0; i<noOptima; i++)
		{
		  sprintf(f + strlen(f),"optimum:%i&nbsp;&nbsp;&nbsp;mean&nbsp;absorb.&nbsp;prob.:%f.2<BR><BR>",optimaIndex[i],meanReachOptFL[i]);
		}
		
		// :: FOURIER DECOMPOSITION BRIEF OUTPUT ::
			sprintf(f + strlen(f),"&nbsp;Fourier Decomposition:<BR>Order\tFraction<BR>");
		if (fourierCoefficient_total!=0)
		{
			for(i=1;i<FL->nlocus+1; i++)
			{
			  sprintf(f + strlen(f),"&nbsp;&nbsp;&nbsp;%d&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;%f<BR>",i,fourierCoefficients.val[i][0]/fourierCoefficient_total);
			}
		}
		else
		{
		  sprintf(f + strlen(f),"none<BR>");
		}
	
	
		if(opt_log)
			exp_landscape( FL );

		free(HistoOutDegree);
		free_chain( mychain );
	
		free(gamma_star);
		free(gamma);
		free(e);
	
		FreeMat(&P);
		FreeMat(&Q);
		FreeMat(&R);
		FreeMat(&N);
		FreeMat(&expectedNoSteps);
		FreeMat(&varianceNoSteps);
		FreeMat(&transientProbs);
		FreeMat(&B);
		FreeMat(&reachibility);
		FreeMat(&fitterGenotypes);
	
		FreeMat(&fourierCoefficients);
	
		FreeMat(&GammaMultiAllele_AiBiAjBj_Denominator);
		FreeMat(&GammaMultiAllele_AiBiAjBj_Nominator);
		FreeMat(&GammaMultiAllele_AiBiAjBj_Matrix);
	
		FreeMat(&GammaMultiAllele_AiBi_Denominator);
		FreeMat(&GammaMultiAllele_AiBi_Nominator);
		FreeMat(&GammaMultiAllele_AiBi_Matrix);
		
		FreeMat(&GammaMultiAllele_AjBj_Denominator);
		FreeMat(&GammaMultiAllele_AjBj_Nominator);
		FreeMat(&GammaMultiAllele_AjBj_Matrix);
	
		FreeMat(&GammaMultiAllele_AiBj_Nominator);
		FreeMat(&GammaMultiAllele_AiBj_Denominator);
		FreeMat(&GammaMultiAllele_AiBj_Matrix);
		
		FreeMat(&GammaMultiAllele_Ai_Denominator);
		FreeMat(&GammaMultiAllele_Ai_Nominator);
		FreeMat(&GammaMultiAllele_Ai_Matrix);
		
		FreeMat(&GammaMultiAllele_Bj_Denominator);
		FreeMat(&GammaMultiAllele_Bj_Nominator);
		FreeMat(&GammaMultiAllele_Bj_Matrix);
	
		FreeMat(&GammaMultiAllele_ij_Denominator);
		FreeMat(&GammaMultiAllele_ij_Nominator);
		FreeMat(&GammaMultiAllele_ij_Matrix);
	
		FreeMat(&GammaMultiAllele_i_Denominator);
		FreeMat(&GammaMultiAllele_i_Nominator);
		FreeMat(&GammaMultiAllele_i_Matrix);

		FreeMat(&GammaMultiAllele_j_Denominator);
		FreeMat(&GammaMultiAllele_j_Nominator);
		FreeMat(&GammaMultiAllele_j_Matrix);
		return f;
	
	}



/*
	outputs some nice values about lanscape
	IncreaseRatio - how much the fitness should be different
	opt_log       - use the log of fitnesses
	opt_short     - can be 0, 1 or 2. 0, expanded output, 1 short report, and 2 skip the first line
	FLname        - a name for the FL
	outfile       - if NULL, output in stdout, otehrwise, append to outfile
	
*/
void OutputSummaryStats( struct landscape * FL, float IncreaseRatio, char opt_log , char opt_short, char *FLname, char *outfile )
{
	int i,j,q,g;
	int noOptima = 0;
	double meanExpectedNumberOfStepsFL = 0;
	double mean2ExpectedNumberOfStepsFL = 0;
	double varExpectedNumberOfStepsFL;
	
	double meanVarianceNumberOfStepsFL = 0;
	double mean2VarianceNumberOfStepsFL = 0;
	double varVarianceNumberOfStepsFL;
	
	double meanReachabilityFL = 0;
	double mean2ReachabilityFL = 0;
	double varReachabilityFL;
	
	double meanFitterGenotypesFL = 0;
	double mean2FitterGenotypesFL = 0;
	// double varFitterGenotyspesFL;
	
	int npeaks;                       /* number of peaks and sinks */
	int nsinks;
	

	/* general measure of FL curvature */
	double * gamma = (double *)malloc( (size_t) FL->nlocus*sizeof(double) );
	double * gamma_star = (double *)malloc( (size_t) FL->nlocus*sizeof(double) );
	
	int neighbors;                     /* how many neighbors for each genotype */
	int *HistoOutDegree=NULL;
	double mean_outdegree, mean2_outdegree, var_outdegree;
	
	struct matrix P = MemMat(FL->ngenotypes, FL->ngenotypes);
	compute_TransitionMatrix(FL, &P, &noOptima);					/* Transition matrix for entire fitness landscape*/
	double varFitterGenotypesFL;
	double meanReachOptFL[noOptima];
	double mean2ReachOptFL[noOptima];
	double varReachOptFL[noOptima];
	for(i = 0; i < noOptima; i++)
	{
		meanReachOptFL[i] = 0;
		mean2ReachOptFL[i] = 0;
	}
	
	struct matrix Q = MemMat(P.r - noOptima, P.c - noOptima);
	compute_QMatrix(&P, &Q);										/* Transition matrix between transitive states (non-optima)*/
	
	int * optimaIndex = (int *)malloc( (size_t) noOptima*sizeof(int) );  //int optimaIndex[noOptima];
	
	struct matrix R = MemMat(P.r - noOptima, noOptima);
	compute_RMatrix(&P, noOptima, optimaIndex, &R);					/* Transition matrix between transitive states -> optima; also returns index of optima */
	
	struct matrix N	= MemMat(P.r - noOptima, P.c - noOptima);
	compute_FundamentalMatrix(&Q, &N);								/* Fundamental matrix N = (I-Q)^-1 */
	
	struct matrix expectedNoSteps = MemMat(N.r, 1);
	compute_ExpectedNumberOfSteps(&N, &expectedNoSteps);			/* Expected number of steps (until being absorbed/ any optimum is reached) t = N 1 */
	
	struct matrix varianceNoSteps = MemMat(N.r, 1);
	compute_VarianceNumberOfSteps(&N, &expectedNoSteps, &varianceNoSteps);	/*Variance in number of steps (until being absorbed/ any optimum is reached) = (2N-I)t - t'*'t ('*': Hadarmard product) */
	
	struct matrix transientProbs = MemMat(N.r, N.c);
	compute_TransientProbs(&N, &transientProbs);					/* Transient probabilities to reach any non-absorbing state from any non-absorbing state (eventually) */
	
	struct matrix B = MemMat(P.r - noOptima, noOptima);
	compute_AbsorbingProbs(&N, &R, &B);								/* Absorbing probabilities to reach any optimum from any non-absorbing state (eventually)*/
	
	struct matrix reachibility = MemMat(Q.r, 1);
	struct matrix fitterGenotypes = MemMat(Q.r, 1);

	i = 0;
	for (g = 0; g < FL->ngenotypes; g++)
	{
		if (P.val[g][g] != 1.)
		{
			reachibility.val[i][0] = computeReachibility(&transientProbs, &B, i);	/* Check if given transient state can be reached from a given starting genotype, i.e., when their transient probability is > 0 */
			fitterGenotypes.val[i][0] = CountFitterGenotypes(FL, g, 1., 1);			/* Count the number of mutants with higher fitness than that of a given starting genotype */
			
			meanExpectedNumberOfStepsFL += expectedNoSteps.val[i][0];
			mean2ExpectedNumberOfStepsFL += expectedNoSteps.val[i][0] * expectedNoSteps.val[i][0];
			
			meanVarianceNumberOfStepsFL += varianceNoSteps.val[i][0];
			mean2VarianceNumberOfStepsFL += varianceNoSteps.val[i][0] * varianceNoSteps.val[i][0];
			
			meanReachabilityFL += reachibility.val[i][0];
			mean2ReachabilityFL += reachibility.val[i][0] * reachibility.val[i][0];
			
			meanFitterGenotypesFL += fitterGenotypes.val[i][0];
			mean2FitterGenotypesFL += fitterGenotypes.val[i][0] *fitterGenotypes.val[i][0];
			
			for (j = 0; j < noOptima; j++)
			{
				meanReachOptFL[j] += B.val[i][j];
				mean2ReachOptFL[j] += B.val[i][j] * B.val[i][j];
			}
			i++;
		}
	}
	
	meanExpectedNumberOfStepsFL /= (0.0 + FL->ngenotypes - noOptima);
	mean2ExpectedNumberOfStepsFL /= (0.0 + FL->ngenotypes - noOptima);
	varExpectedNumberOfStepsFL = ((FL->ngenotypes - noOptima)/((FL->ngenotypes - noOptima)-1.0))*( mean2ExpectedNumberOfStepsFL - pow(meanExpectedNumberOfStepsFL,2.0) );
	
	meanVarianceNumberOfStepsFL /= (0.0 + FL->ngenotypes - noOptima);
	mean2VarianceNumberOfStepsFL /= (0.0 + FL->ngenotypes - noOptima);
	varVarianceNumberOfStepsFL = ((FL->ngenotypes - noOptima)/((FL->ngenotypes - noOptima)-1.0))*( mean2VarianceNumberOfStepsFL - pow(meanVarianceNumberOfStepsFL,2.0) );
	
	meanReachabilityFL /= (0.0 + FL->ngenotypes - noOptima);
	mean2ReachabilityFL /= (0.0 + FL->ngenotypes - noOptima);
	varReachabilityFL = ((FL->ngenotypes - noOptima)/((FL->ngenotypes - noOptima)-1.0))*( mean2ReachabilityFL - pow(meanReachabilityFL,2.0) );
	
	meanFitterGenotypesFL /= (0.0 + FL->ngenotypes - noOptima);
	mean2FitterGenotypesFL /= (0.0 + FL->ngenotypes - noOptima);
	
	varFitterGenotypesFL = ((FL->ngenotypes - noOptima)/((FL->ngenotypes - noOptima)-1.0))*( mean2FitterGenotypesFL - pow(meanFitterGenotypesFL,2.0) );
	
	for (j = 0; j < noOptima; j++)
	{
		meanReachOptFL[j] /= (0.0 + FL->ngenotypes - noOptima);
		mean2ReachOptFL[j] /= (0.0 + FL->ngenotypes - noOptima);
		varReachOptFL[j] = ((FL->ngenotypes - noOptima)/((FL->ngenotypes - noOptima)-1.0))*( mean2ReachOptFL[j] - pow(meanReachOptFL[j],2.0) );
	}
	
	
	int sizeSubSets = 0;
	int noAlleles = 0;
	for (j = 0; j < FL->nlocus; j++)
	{
		sizeSubSets += N_choose_n(FL->alleles[j], 2);
		noAlleles += FL->alleles[j];
	}
	
	
	struct matrix GammaMultiAllele_AiBiAjBj_Nominator = MemMat(sizeSubSets, sizeSubSets);
	struct matrix GammaMultiAllele_AiBiAjBj_Denominator = MemMat(sizeSubSets, sizeSubSets);
	struct matrix GammaMultiAllele_AiBiAjBj_Matrix = MemMat(sizeSubSets, sizeSubSets);
	
	struct matrix GammaMultiAllele_AiBi_Nominator = MemMat(sizeSubSets, 1);
	struct matrix GammaMultiAllele_AiBi_Denominator = MemMat(sizeSubSets, 1);
	struct matrix GammaMultiAllele_AiBi_Matrix = MemMat(sizeSubSets, 1);
	
	struct matrix GammaMultiAllele_AiBj_Nominator = MemMat(noAlleles, noAlleles);
	struct matrix GammaMultiAllele_AiBj_Denominator = MemMat(noAlleles, noAlleles);
	struct matrix GammaMultiAllele_AiBj_Matrix = MemMat(noAlleles, noAlleles);
	
	struct matrix GammaMultiAllele_AjBj_Nominator = MemMat(sizeSubSets, 1);
	struct matrix GammaMultiAllele_AjBj_Denominator = MemMat(sizeSubSets, 1);
	struct matrix GammaMultiAllele_AjBj_Matrix = MemMat(sizeSubSets, 1);
	
	struct matrix GammaMultiAllele_ij_Nominator = MemMat(FL->nlocus, FL->nlocus);
	struct matrix GammaMultiAllele_ij_Denominator = MemMat(FL->nlocus, FL->nlocus);
	struct matrix GammaMultiAllele_ij_Matrix = MemMat(FL->nlocus, FL->nlocus);

	
	GammaMultiAllele_AiBiAjBj( *FL, &GammaMultiAllele_AiBiAjBj_Nominator, &GammaMultiAllele_AiBiAjBj_Denominator, &GammaMultiAllele_AiBi_Nominator, &GammaMultiAllele_AiBi_Denominator, &GammaMultiAllele_AjBj_Nominator, &GammaMultiAllele_AjBj_Denominator, &GammaMultiAllele_ij_Nominator, &GammaMultiAllele_ij_Denominator,-1);
	for (i = 0; i < sizeSubSets; i++)
	{
		for (j = 0; j < sizeSubSets; j++)
		{
			GammaMultiAllele_AiBiAjBj_Matrix.val[i][j] = GammaMultiAllele_AiBiAjBj_Nominator.val[i][j] / GammaMultiAllele_AiBiAjBj_Denominator.val[i][j];
		}
		GammaMultiAllele_AiBi_Matrix.val[i][0] = GammaMultiAllele_AiBi_Nominator.val[i][0] / GammaMultiAllele_AiBi_Denominator.val[i][0];
		GammaMultiAllele_AjBj_Matrix.val[i][0] = GammaMultiAllele_AjBj_Nominator.val[i][0] / GammaMultiAllele_AjBj_Denominator.val[i][0];
	}
	
	struct matrix GammaMultiAllele_Ai_Nominator = MemMat(noAlleles, 1);
	struct matrix GammaMultiAllele_Ai_Denominator = MemMat(noAlleles, 1);
	struct matrix GammaMultiAllele_Ai_Matrix = MemMat(noAlleles, 1);
	
	struct matrix GammaMultiAllele_Bj_Nominator = MemMat(noAlleles, 1);
	struct matrix GammaMultiAllele_Bj_Denominator = MemMat(noAlleles, 1);
	struct matrix GammaMultiAllele_Bj_Matrix = MemMat(noAlleles, 1);
	
	GammaMultiAllele_Ai_Bj( *FL, noAlleles, &GammaMultiAllele_AiBi_Nominator, &GammaMultiAllele_AiBi_Denominator, &GammaMultiAllele_AjBj_Nominator, &GammaMultiAllele_AjBj_Denominator, &GammaMultiAllele_Ai_Nominator, &GammaMultiAllele_Ai_Denominator, &GammaMultiAllele_Bj_Nominator, &GammaMultiAllele_Bj_Denominator, -1);
	GammaMultiAllele_AiBj( *FL, noAlleles, &GammaMultiAllele_AiBiAjBj_Nominator, &GammaMultiAllele_AiBiAjBj_Denominator, &GammaMultiAllele_AiBj_Nominator, &GammaMultiAllele_AiBj_Denominator, -1);
	for (i = 0; i < noAlleles; i++)
	{
		GammaMultiAllele_Ai_Matrix.val[i][0] = GammaMultiAllele_Ai_Nominator.val[i][0] / GammaMultiAllele_Ai_Denominator.val[i][0];
		GammaMultiAllele_Bj_Matrix.val[i][0] = GammaMultiAllele_Bj_Nominator.val[i][0] / GammaMultiAllele_Bj_Denominator.val[i][0];
		for (j = 0; j < noAlleles; j++)
		{
			GammaMultiAllele_AiBj_Matrix.val[i][j] = GammaMultiAllele_AiBj_Nominator.val[i][j] / GammaMultiAllele_AiBj_Denominator.val[i][j];
		}
	}
	
	struct matrix GammaMultiAllele_i_Nominator = MemMat(FL->nlocus, 1);
	struct matrix GammaMultiAllele_i_Denominator = MemMat(FL->nlocus, 1);
	struct matrix GammaMultiAllele_i_Matrix = MemMat(FL->nlocus, 1);
	
	struct matrix GammaMultiAllele_j_Nominator = MemMat(FL->nlocus, 1);
	struct matrix GammaMultiAllele_j_Denominator = MemMat(FL->nlocus, 1);
	struct matrix GammaMultiAllele_j_Matrix = MemMat(FL->nlocus, 1);
	
	
	GammaMultiAllele_ij( *FL, noAlleles, &GammaMultiAllele_Ai_Nominator, &GammaMultiAllele_Ai_Denominator, &GammaMultiAllele_Bj_Nominator, &GammaMultiAllele_Bj_Denominator, &GammaMultiAllele_i_Nominator, &GammaMultiAllele_i_Denominator, &GammaMultiAllele_j_Nominator, &GammaMultiAllele_j_Denominator, -1);
	for (i = 0; i < FL->nlocus; i++)
	{
		for (j = 0; j < FL->nlocus; j++)
		{
			GammaMultiAllele_ij_Matrix.val[i][j] = GammaMultiAllele_ij_Nominator.val[i][j] / GammaMultiAllele_ij_Denominator.val[i][j];
		}
		GammaMultiAllele_i_Matrix.val[i][0] = GammaMultiAllele_i_Nominator.val[i][0] / GammaMultiAllele_i_Denominator.val[i][0];
		GammaMultiAllele_j_Matrix.val[i][0] = GammaMultiAllele_j_Nominator.val[i][0] / GammaMultiAllele_j_Denominator.val[i][0];
	}
	
	float rs;

	float *HistoPairs;
	int npairs;

	double *moments=NULL;
	float *array=NULL;

	FILE *f;
	
	double fourierCoefficient_total = 0;
	double fourierCoefficient_epi = 0;
	
	struct matrix fourierCoefficients = MemMat(FL->nlocus + 1, 1);
	int max_f = 1;
	
	for (i = 0; i < FL->nlocus + 1; i++)
	{
		fourierCoefficients.val[i][0] = 0;
	}
	FourierTransform(FL, &fourierCoefficients, &fourierCoefficient_total, &fourierCoefficient_epi, &max_f);
	
	int *e, total_e;
	
	
	struct Chains mychain;            /* structure to store chain results */
	long nsteps=0,
	    norigins=0,
	    depthmax=0;

	
	if(opt_log)
	{
		log_landscape( FL );
	}
	
	if(outfile)
	{
		f=fopen(outfile, "a");
		if(!f)fprintf(stderr, "cannot write into file %s, please check. bye\n", outfile), exit(3);
	}
	else
	{
		f=stdout;
	}
	
	neighbors = 0;
	for( i=0; i<FL->nlocus; i++)
	{
		neighbors += FL->alleles[i] - 1;
	}

	
	npeaks = numberPeaks( FL, IncreaseRatio );
	nsinks = numberSinks( FL, IncreaseRatio );
	gamma[0] = 1;
	gamma_star[0] = 1;
	for (i = 1; i < FL->nlocus; i++)
	{
		gamma[i] = GammaDistance(*FL, i, -1);
		gamma_star[i] = GammaDistance(*FL, i, 0);
	}
	
	HistoOutDegree = compute_HistoBetterNeighbor( FL, IncreaseRatio );
	mean_outdegree=0;
	mean2_outdegree=0;

	for(i=0;i<=neighbors;i++)
	{
		mean_outdegree += i*HistoOutDegree[i];
		mean2_outdegree += i*i*HistoOutDegree[i];
	}
	mean_outdegree/=(0.0+FL->ngenotypes);
	mean2_outdegree/=(0.0+FL->ngenotypes);
	var_outdegree=(FL->ngenotypes/(FL->ngenotypes-1.0))*( mean2_outdegree - pow(mean_outdegree,2.0) );
	
	mychain = Compute_Chains( FL, 0, IncreaseRatio  );
	
	rs = compute_rs(  FL );
	
	e = get_sign_epistasis_FL( FL );
	total_e=e[0]+e[1]+e[2]+e[3];

	for( i=0 ; i < mychain.nchains ; i++ ){
		nsteps += mychain.steps[i];
		norigins += mychain.origins[i];
		depthmax = (mychain.depth[i]<depthmax)?depthmax:mychain.depth[i];
	}
	
	HistoPairs = Count_AllLocusPairs( FL, IncreaseRatio  );
	npairs = FL->nlocus*(FL->nlocus-1);
	
	
	array = (float *)malloc( (size_t) npairs*sizeof(float) );
	if(!array)fprintf(stderr, "cannot allocated array, bye\n"), exit(1);
	q = 0;
	
	for(i=0;i<FL->nlocus;i++)
	{
		for(j=0;j<FL->nlocus;j++)
		{
			if(i!=j)
			{
				array[q++] = HistoPairs[ i*FL->nlocus + j ] ; /*  +HistoPairs[ j*h.nlocus + i ];*/
			}
		}
	}
	
	moments = Moments( array , npairs, 'f' );
	
	

	/****************************
		::   Brief Output   ::
	 ****************************/
	if(opt_short)
	{
		for(i=0;i<23+noOptima;i++)
		{
			fprintf(f,"(%d)\t",i+1);
		}
		fprintf(f,"\n");
		fprintf(f,"name\tngeno\tnpeaks\tnsinks\tgamma\tgamma*\tr/s\tnchains\tnsteps\tnori\tdepth\t"
				"magn\tsign\trsign\tf[1]\t[2]\tf[3+]\tmode_f\toutD_m\toutD_v\tsteps_m\treach_m\tfitG_m");
		
		for (i=0; i<noOptima; i++)
		{
			fprintf(f,"\topt_i\tmProbOpt_%i",i);
		}
		fprintf(f,"\n");
		
		if (strchr(FLname,'/') == 0 )
			fprintf(f,"%s\t%d\t%d\t%d\t", FLname, FL->ngenotypes, npeaks, nsinks);
		else
		{		
			int tmp=  ( strchr(FLname,'/') - FLname );
			FLname[tmp]=0;
			fprintf(f,"%s\t%d\t%d\t%d\t", FLname, FL->ngenotypes, npeaks, nsinks);
			FLname[tmp]='/';

		}

		
		fprintf(f,"%.3f\t%.3f\t%.3f\t", gamma[1], gamma_star[1], rs);
		fprintf(f,"%ld\t%ld\t%ld\t%ld\t",  mychain.nchains, nsteps, norigins, depthmax);
		fprintf(f,"%.3f\t%.3f\t%.3f\t", e[1]/(0.0+total_e), e[2]/(0.0+total_e), e[3]/(0.0+total_e) );
		
		fprintf(f,"%.3f\t%.3f\t%.3f\t%d\t", fourierCoefficients.val[1][0] / fourierCoefficient_total,  fourierCoefficients.val[2][0] / fourierCoefficient_total, 1 - ( fourierCoefficients.val[1][0] + fourierCoefficients.val[2][0] )/fourierCoefficient_total, max_f );
		
		fprintf(f,"%.3f\t%.3f\t", sqrt(moments[1]), var_outdegree );
		fprintf(f,"%.3f\t%.3f\t%.3f\t", meanExpectedNumberOfStepsFL, meanReachabilityFL, meanFitterGenotypesFL );

		
		
		for (i = 1; i < FL->nlocus + 1; i++)
		{
			printf("F[%d] = %.3f\n", i, fourierCoefficients.val[i][0] / fourierCoefficient_total);
		}

		
	
		for (i=0; i<noOptima; i++)
		{
			fprintf(f,"%d\t%.3f\t",optimaIndex[i], meanReachOptFL[i]);
		}
		fprintf(f,"\n");
		
		
		if(opt_log)
		{
			exp_landscape( FL );
		}
		
		free(HistoOutDegree);
		free_chain( mychain );
		
		free(gamma_star);
		free(gamma);
		free(e);
		
		free(optimaIndex);
		FreeMat(&P);
		FreeMat(&Q);
		FreeMat(&R);
		FreeMat(&N);
		FreeMat(&expectedNoSteps);
		FreeMat(&varianceNoSteps);
		FreeMat(&transientProbs);
		FreeMat(&B);
		FreeMat(&reachibility);
		FreeMat(&fitterGenotypes);
		
		FreeMat(&fourierCoefficients);
		
		FreeMat(&GammaMultiAllele_AiBiAjBj_Denominator);
		FreeMat(&GammaMultiAllele_AiBiAjBj_Nominator);
		FreeMat(&GammaMultiAllele_AiBiAjBj_Matrix);
		
		FreeMat(&GammaMultiAllele_AiBi_Denominator);
		FreeMat(&GammaMultiAllele_AiBi_Nominator);
		FreeMat(&GammaMultiAllele_AiBi_Matrix);
		
		FreeMat(&GammaMultiAllele_AiBj_Nominator);
		FreeMat(&GammaMultiAllele_AiBj_Denominator);
		FreeMat(&GammaMultiAllele_AiBj_Matrix);
		
		FreeMat(&GammaMultiAllele_AjBj_Denominator);
		FreeMat(&GammaMultiAllele_AjBj_Nominator);
		FreeMat(&GammaMultiAllele_AjBj_Matrix);
		
		FreeMat(&GammaMultiAllele_Ai_Denominator);
		FreeMat(&GammaMultiAllele_Ai_Nominator);
		FreeMat(&GammaMultiAllele_Ai_Matrix);
		
		FreeMat(&GammaMultiAllele_Bj_Denominator);
		FreeMat(&GammaMultiAllele_Bj_Nominator);
		FreeMat(&GammaMultiAllele_Bj_Matrix);
		
		FreeMat(&GammaMultiAllele_ij_Denominator);
		FreeMat(&GammaMultiAllele_ij_Nominator);
		FreeMat(&GammaMultiAllele_ij_Matrix);
		
		FreeMat(&GammaMultiAllele_i_Denominator);
		FreeMat(&GammaMultiAllele_i_Nominator);
		FreeMat(&GammaMultiAllele_i_Matrix);
		
		FreeMat(&GammaMultiAllele_j_Denominator);
		FreeMat(&GammaMultiAllele_j_Nominator);
		FreeMat(&GammaMultiAllele_j_Matrix);
		
		return;
	
	}
	

	/***************************
		:: Expanded Output ::
	****************************/

	fprintf(f,"/* FL name */\n");
	fprintf(f,"   %s\n", FLname);


	/*
		Print out Peaks/Sinks
	*/
	fprintf(f,"\n/* Peaks/Sinks */\n");
	
	fprintf(f,"   #genotypes: %d\n", FL->ngenotypes );
	fprintf(f,"   #peaks: %d\n", npeaks);
	fprintf(f,"   #sinks: %d\n", nsinks);
	
	/*
		Print out the fraction of pairwise epistasis
	*/
	fprintf(f,"\n/* Epistasis types */\n");
	
	fprintf(f,"   none:       %.3f\n", e[0]/(0.0+total_e) );
	fprintf(f,"   magnitude:  %.3f\n", e[1]/(0.0+total_e));
	fprintf(f,"   sign:       %.3f\n", e[2]/(0.0+total_e));
	fprintf(f,"   reciprocal: %.3f\n", e[3]/(0.0+total_e));
	
	
	
	/*
		Roughness / Slope
	*/
	
	fprintf(f,"\n/* Roughness / Slope */\n");
	fprintf(f,"   r/s: %f\n", rs  );

	/*
		Gamma
	*/
	
	fprintf(f,"\n\n");
	
	fprintf(f,"/*****************/\n");
	fprintf(f,"/****  Gamma  ****/\n");
	fprintf(f,"/*****************/\n");

	fprintf(f,"\n/* Gamma Global */\n");

	fprintf(f,"   gamma[0]:\t1\n");
	for (i = 1; i < FL->nlocus; i++)
	{
		fprintf(f,"   gamma[%d]:\t%f\n", i, gamma[i]);
	}
	fprintf(f,"\n");

	fprintf(f,"\n/* Gamma Locus */\n\n");
	for (i = 0; i < FL->nlocus; i++)
	{
		fprintf(f,"   Gamma(%d->):\t%f\n", i+1, GammaMultiAllele_j_Matrix.val[i][0]);
	}
	fprintf(f,"\n");
	
	for (i = 0; i < FL->nlocus; i++)
	{
		fprintf(f,"   Gamma(->%d):\t%f\n", i+1, GammaMultiAllele_i_Matrix.val[i][0]);
	}
	fprintf(f,"\n");
	

	for (i = 0; i < FL->nlocus; i++)
	{
		for (j = 0; j < FL->nlocus; j++)
		{
			fprintf(f,"   Gamma(%d->%d):\t%f\n", i+1, j+1, GammaMultiAllele_ij_Matrix.val[j][i]);
		}
		fprintf(f,"\n");
	}
	fprintf(f,"\n");
	
	
	fprintf(f,"\n/* Gamma Allele Pairs */\n\n");
	q = 1;
	for (i = 0; i < sizeSubSets; i++)
	{
		g = 1;
		for (j = 0; j < sizeSubSets; j++)
		{
			fprintf(f,"   Gamma(%d->%d):\t%f\n", q, g, GammaMultiAllele_AiBiAjBj_Matrix.val[j][i]);
			g++;
		}
		q++;
		fprintf(f,"\n");
	}
	fprintf(f,"\n");

	for (j = 0; j < sizeSubSets; j++)
	{
		fprintf(f,"   Gamma(->%d):\t%f\n", j+1, GammaMultiAllele_AiBi_Matrix.val[j][0]);
	}
	fprintf(f,"\n");

	
	for (j = 0; j < sizeSubSets; j++)
	{
		fprintf(f,"   Gamma(%d->):\t%f\n", j+1, GammaMultiAllele_AjBj_Matrix.val[j][0]);
	}
	fprintf(f,"\n");
	
	fprintf(f,"\n/* Gamma Allele */\n\n");

	for (i = 0; i < noAlleles; i++)
	{
		for (j = 0; j < noAlleles; j++)
		{
			fprintf(f,"   Gamma(%d->%d):\t%f\n", i+1, j+1, GammaMultiAllele_AiBj_Matrix.val[j][i]);
		}
		fprintf(f,"\n");
	}
	fprintf(f,"\n");

	
	for (j = 0; j < noAlleles; j++)
	{
		fprintf(f,"   Gamma(->%d):\t%f\n", j+1, GammaMultiAllele_Ai_Matrix.val[j][0]);
	}
	fprintf(f,"\n");

	for (j = 0; j < noAlleles; j++)
	{
		fprintf(f,"   Gamma(%d->):\t%f\n", j+1, GammaMultiAllele_Bj_Matrix.val[j][0]);
	}
	fprintf(f,"\n");


	/*
		FOURIER DECOMPOSITION
	*/

	
	fprintf(f,"\n/* Fourier Decomposition */\n");
	fprintf(f,"   o\tfraction\n" );
	
	for(i=1;i<FL->nlocus+1; i++)
	{

			fprintf(f,"   %d\t%.3f\n", i, fourierCoefficients.val[i][0] / fourierCoefficient_total);
	}
	
	fprintf(f,"\n/* Fourier Fraction epistasis */\n" );
	fprintf(f,"\t%.3f\n", fourierCoefficient_epi/fourierCoefficient_total);
	
	
	/*
		Print out steps
	*/


	fprintf(f,"\n/* Out-degree distribution */\n");
	
	fprintf(f,"   deg #   fract\n");
	for(i=0;i<=neighbors;i++){
		fprintf(f, "d:  %d  %3d %.3f\n", i, HistoOutDegree[ i ], HistoOutDegree[ i ] / (float) FL->ngenotypes );

	}
	
	fprintf(f,"Mean(OutDegree):  %f\n", mean_outdegree );
	fprintf(f,"Stdev(OutDegree): %f\n", var_outdegree );


	/*
		see what genotype graph has to say
	*/


	fprintf(f,"\n/* Accessible genotypes from g after (s steps) */\n");
	fprintf(f,"   g\ts1   s2   s3 ...\n");
	fprintf(f,"----\t----------------\n");

	int *tmpg=NULL;
	for(q=0;q<FL->ngenotypes;q++)
	{
	
		struct list res;
		
		res = count_genotypes( q, FL, 0, 0, 0, IncreaseRatio );   /* count how many genotypes are reached at distance of X, assuming constant increase */
		
		tmpg = int2genotype( *FL, q , tmpg);
		print_genotype( tmpg, FL->nlocus );
		fprintf(f,"   %d\t 1 ",q);
		
		
		for(i=0; i<res.n; i++){
			fprintf(f,"%4ld ", res.genotypes[i] );
		}
		fprintf(f,"\n");
		
	}
	
	free(tmpg);
	

	/*
	 
	   CHAINS
	 
	*/

	fprintf(f,"\n/* Chains */\n");
	fprintf(f," nchains: %ld\n",  mychain.nchains);

	fprintf(f,"   id\tsteps\tdepth\torigins\n");
	for( i=0 ; i < mychain.nchains ; i++ )	
		fprintf(f,"   %d\t%ld\t%ld\t%ld\n",  i+1, mychain.steps[i],  mychain.depth[i],  mychain.origins[i]);
	
	free_chain( mychain );
	
	/*
	
		Pairs of Loci
	 
	*/

	fprintf(f,"\n/* Count when mu at loci i->j are both beneficials (first i, then j) */\n");

	q = 0;
	for(i=0;i<FL->nlocus;i++)
		for(j=0;j<FL->nlocus;j++)
			if(i!=j)
			{
				fprintf(f, "   %d->%d: %d\n", i+1, j+1, (int) array[q++] );
				
			}
			
	fprintf(f,"Pairs Mean: %f ; Stdev: %f ; Skewness: %f\n", moments[0], sqrt(moments[1]), moments[2] );


	/*
	 
		ADAPTIVE WALKS
	 
	 */
	
	fprintf(f,"\n/* Adaptive walk statistics */\n");
	fprintf(f,"G\tE[steps]\tV[steps]\tReachblty\tFitterG");
	for (i = 0; i < noOptima; i++)
	{
		fprintf(f,"\tP[%d]",optimaIndex[i]);
	}
	fprintf(f,"\n");

	q=0;
	for (i = 0; i < FL->ngenotypes; i++)
	{
		if(P.val[i][i] != 1)
		{
			fprintf(f, "%d\t%.3f\t%.3f\t%.0f\t%.0f\t", i, expectedNoSteps.val[q][0], varianceNoSteps.val[q][0], reachibility.val[q][0], fitterGenotypes.val[q][0]);
		
			for (j = 0; j < noOptima; j++)
			{
				fprintf(f,"\t%.3f",B.val[q][j]);
			}
			fprintf(f,"\n");
			q++;
		}
	}
	
	fprintf(f,"\n");
	fprintf(f,"\n/* Summary adaptive walk statistics */\n");
	fprintf(f,"Mean E[steps]: %.3f\tStdev E[steps]: %.3f\tMean V[steps]: %.3f\tStdev V[steps]: %.3f\tMean Reachability: %.3f\tStdev Reachability: %.3f\tMean FitterG: %.3f\tStdev FitterG: %.3f",meanExpectedNumberOfStepsFL, sqrt(varExpectedNumberOfStepsFL), meanVarianceNumberOfStepsFL, sqrt(varVarianceNumberOfStepsFL), meanReachabilityFL, sqrt(varReachabilityFL), meanFitterGenotypesFL, sqrt(varFitterGenotypesFL));
	for (i = 0; i < noOptima; i++)
	{
		fprintf(f,"\tMean P[%d]: %.3f\tVar P[%d]: %.3f",optimaIndex[i], meanReachOptFL[i], optimaIndex[i], sqrt(varReachOptFL[i]));
	}
	fprintf(f,"\n");
	
	
	if(opt_log)
	{
		exp_landscape( FL );
	}

	if(outfile)
	{
		fclose(f);
	}
	
	free( HistoPairs );


	free(array);
	free(moments);
	free(optimaIndex);
	free(e);
	free(gamma_star);
	free(gamma);
	
	
	FreeMat(&P);
	FreeMat(&Q);
	FreeMat(&R);
	FreeMat(&N);
	FreeMat(&expectedNoSteps);
	FreeMat(&varianceNoSteps);
	FreeMat(&transientProbs);
	FreeMat(&B);
	FreeMat(&reachibility);
	FreeMat(&fitterGenotypes);
	
	FreeMat(&fourierCoefficients);
	
	FreeMat(&GammaMultiAllele_AiBiAjBj_Denominator);
	FreeMat(&GammaMultiAllele_AiBiAjBj_Nominator);
	FreeMat(&GammaMultiAllele_AiBiAjBj_Matrix);
	
	FreeMat(&GammaMultiAllele_AiBi_Denominator);
	FreeMat(&GammaMultiAllele_AiBi_Nominator);
	FreeMat(&GammaMultiAllele_AiBi_Matrix);
	
	FreeMat(&GammaMultiAllele_AjBj_Denominator);
	FreeMat(&GammaMultiAllele_AjBj_Nominator);
	FreeMat(&GammaMultiAllele_AjBj_Matrix);
	
	FreeMat(&GammaMultiAllele_AiBj_Nominator);
	FreeMat(&GammaMultiAllele_AiBj_Denominator);
	FreeMat(&GammaMultiAllele_AiBj_Matrix);
	
	FreeMat(&GammaMultiAllele_Ai_Denominator);
	FreeMat(&GammaMultiAllele_Ai_Nominator);
	FreeMat(&GammaMultiAllele_Ai_Matrix);
	
	FreeMat(&GammaMultiAllele_Bj_Denominator);
	FreeMat(&GammaMultiAllele_Bj_Nominator);
	FreeMat(&GammaMultiAllele_Bj_Matrix);
	
	FreeMat(&GammaMultiAllele_ij_Denominator);
	FreeMat(&GammaMultiAllele_ij_Nominator);
	FreeMat(&GammaMultiAllele_ij_Matrix);
	
	FreeMat(&GammaMultiAllele_i_Denominator);
	FreeMat(&GammaMultiAllele_i_Nominator);
	FreeMat(&GammaMultiAllele_i_Matrix);
	
	FreeMat(&GammaMultiAllele_j_Denominator);
	FreeMat(&GammaMultiAllele_j_Nominator);
	FreeMat(&GammaMultiAllele_j_Matrix);

	free( HistoOutDegree );
	
}


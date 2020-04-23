/*
	correlation functions at distance d = 1
	only binary landscapes
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "landscape.h"
#include "genotypes.h"
#include "summary_statistics.h"
#include "calculus.h"



/*
	Swicthes the allele at position i : 0 <-> 1
*/
int gmod(struct landscape h, int g, int i)
{
	
	/* change byte i-1 of g */
	
	int l = h.nlocus;
	int genotype[l];
	
	int2genotype(h,g,genotype);
	genotype[i]=1-genotype[i];
	
	return (genotype2int(h,genotype));
}

int gmodMultiAllele(struct landscape h, int g, int locus, int allele)
{
	
	/* change byte i of g */
	
	int l = h.nlocus;
	int genotype[l];
	
	int2genotype(h,g,genotype);
	genotype[locus]=allele;
	
	return (genotype2int(h,genotype));
}


/*
	Return the fitness difference between genotype g and genotype g with
	the other allele at locus j.

	if tolerance is >0, it computes gamma * with some tolerance
*/
double deltaf(struct landscape h, int g, int j, double tolerance)
{

	if (tolerance<0)
	{
		return (h.fitness[gmod(h,g,j)]-h.fitness[g]);
	}
	else 
	{
		double df=h.fitness[gmod(h,g,j)]-h.fitness[g];
		
		if (df>tolerance)
		{
			return (1.0);
		}
		else
		{
			if (df>-tolerance)
			{
				return (0.0);
			}
			else
			{
				return (-1.0);
			}
		}
	}
}

/*
	Return the fitness difference between genotype g and genotype g with
	another 'allele' at locus j.
 
	if tolerance is >0, it computes gamma * with some tolerance
 */

double deltafMultiAllele(struct landscape h, int g, int locus, int allele, double tolerance)
{
	
	if (tolerance<0)
	{
		return (h.fitness[gmodMultiAllele(h, g, locus, allele)]-h.fitness[g]); //
	}
	else
	{
		double df=h.fitness[gmodMultiAllele(h, g, locus, allele)]-h.fitness[g];
		
		if (df>tolerance)
		{
			return (1.0);
		}
		else
		{
			if (df>-tolerance)
			{
				return (0.0);
			}
			else
			{
				return (-1.0);
			}
		}
	}
}

/*
	Caclulates the gamma matrix between pairs of alleles AiBi -> AjBj.
	if tolerance is >0, it computes gamma * with some tolerance
 
*/


void GammaMultiAllele_AiBiAjBj(struct landscape fl, struct matrix * NominatorAiBiAjBj, struct matrix * DenominatorAiBiAjBj, struct matrix * NominatorAiBi, struct matrix * DenominatorAiBi, struct matrix * NominatorAjBj, struct matrix * DenominatorAjBj, struct matrix * Nominatorij, struct matrix * Denominatorij , double tolerance)
{
	int a1,a2,b1,b2,g,lA,lB;
	double s_j = 0, s_ij = 0;
	int *refGeno = NULL;
	int mutI = 0;
	
	int countIndex1 = 0, countIndex2;
	
	//initialize
	for(lA = 0; lA < fl.nlocus; lA++)
	{
		for(lB = 0; lB < fl.nlocus; lB++)
		{
			if (lA ==lB)
			{
				Nominatorij->val[lA][lB] = 1.;
				Denominatorij->val[lA][lB] = 1.;
			}
			else
			{
				Nominatorij->val[lA][lB] = 0.;
				Denominatorij->val[lA][lB] = 0.;
			}
		}
	}
	
	for(lA = 0; lA < fl.nlocus; lA++)
	{
		//loop through all subsets of alleles of size 2 at locus lA
		for(a1 = 0; a1 < fl.alleles[lA]; a1++)
		{
			for(a2 = a1+1; a2 < fl.alleles[lA]; a2++)
			{
				countIndex2 = 0;
				NominatorAiBi->val[countIndex1][0] = 0.;
				DenominatorAiBi->val[countIndex1][0] = 0.;
				for(lB = 0; lB < fl.nlocus; lB++)
				{
					//loop through all subsets of alleles of size 2 at locus lB
					for(b1 = 0; b1 < fl.alleles[lB]; b1++)
					{
						for(b2 = b1+1; b2 < fl.alleles[lB]; b2++)
						{
							if(countIndex1 == 0)
							{
								NominatorAjBj->val[countIndex2][0] = 0.;
								DenominatorAjBj->val[countIndex2][0] = 0.;
							}
							
							//per Definition gamma = 1 for interactions at the same locus
							if(lA == lB)
							{
								NominatorAiBiAjBj->val[countIndex1][countIndex2] = 1.;
								DenominatorAiBiAjBj->val[countIndex1][countIndex2] = 1.;
								
								Nominatorij->val[lA][lB] = 1.;
								Denominatorij->val[lA][lB] = 1.;

							}
							else
							{
								NominatorAiBiAjBj->val[countIndex1][countIndex2] = 0.;
								DenominatorAiBiAjBj->val[countIndex1][countIndex2] = 0.;

								//loop over all genotypes ...
								for(g = 0; g < fl.ngenotypes; g++)
								{
									refGeno = int2genotype(fl, g, refGeno);
									
									//... but only consider those carrying one of the two focal alleles at locus lA *AND* lB, respectively
									if( ((refGeno[lA] == a1) || (refGeno[lA] == a2)) && ((refGeno[lB] == b1) || (refGeno[lB] == b2)) )
									{
										//create mutant genotypes and calculate corresponding selection coefficient
										if(refGeno[lA] == a1)
										{
											s_j = deltafMultiAllele(fl, g, lA, a2, tolerance);

											if(refGeno[lB] == b1)
											{
												mutI = gmodMultiAllele(fl, g, lB, b2);
												s_ij = deltafMultiAllele(fl, mutI, lA, a2, tolerance);
											}
											else
											{
												mutI = gmodMultiAllele(fl, g, lB, b1);
												s_ij = deltafMultiAllele(fl, mutI, lA, a2, tolerance);
											}
											
										}
										else
										{
											s_j = deltafMultiAllele(fl, g, lA, a1, tolerance);
											
											if(refGeno[lB] == b1)
											{
												mutI = gmodMultiAllele(fl, g, lB, b2);
												s_ij = deltafMultiAllele(fl, mutI, lA, a1, tolerance);
											}
											else
											{
												mutI = gmodMultiAllele(fl, g, lB, b1);
												s_ij = deltafMultiAllele(fl, mutI, lA, a1, tolerance);
											}
											
										}
										
										NominatorAiBiAjBj->val[countIndex1][countIndex2] += s_j * s_ij;
										DenominatorAiBiAjBj->val[countIndex1][countIndex2] += pow(s_j, 2.);
									}
									
								}
								NominatorAiBi->val[countIndex1][0] += NominatorAiBiAjBj->val[countIndex1][countIndex2];
								DenominatorAiBi->val[countIndex1][0] += DenominatorAiBiAjBj->val[countIndex1][countIndex2];
								
								NominatorAjBj->val[countIndex2][0] += NominatorAiBiAjBj->val[countIndex1][countIndex2];
								DenominatorAjBj->val[countIndex2][0] += DenominatorAiBiAjBj->val[countIndex1][countIndex2];
								
								Nominatorij->val[lA][lB] += NominatorAiBiAjBj->val[countIndex1][countIndex2];
								Denominatorij->val[lA][lB] += DenominatorAiBiAjBj->val[countIndex1][countIndex2];

							}
							countIndex2++;
						}
					}
				}
				countIndex1++;
			}
		}
	}
	free(refGeno);
}

/*
	Caclulates the gamma matrix for a focal allele Ai and Bj (i.e., gamma of a mutation involving allele Ai on other mutations and gamma of other mutations on mutations involving Bj).
	if tolerance is >0, it computes gamma * with some tolerance
 
 */

void GammaMultiAllele_Ai_Bj(struct landscape fl, int noAlleles, struct matrix * NominatorAiBi, struct matrix * DenominatorAiBi, struct matrix * NominatorAjBj, struct matrix * DenominatorAjBj, struct matrix * NominatorAi, struct matrix * DenominatorAi, struct matrix * NominatorBj, struct matrix * DenominatorBj, double tolerance)
{
	int a1,a2,lA;
	int alleleCount = 0, indexCount = 0;

	// initialize
	for(a1 = 0; a1 < noAlleles; a1++)
	{
		NominatorAi->val[a1][0] = 0.;
		DenominatorAi->val[a1][0] = 0.;
		NominatorBj->val[a1][0] = 0.;
		DenominatorBj->val[a1][0] = 0.;
	}
	
	for(lA = 0; lA < fl.nlocus; lA++)
	{
		for(a1 = 0; a1 < fl.alleles[lA]; a1++)
		{
			for(a2 = a1+1; a2 < fl.alleles[lA]; a2++)
			{
				NominatorAi->val[alleleCount + a1][0] += NominatorAiBi->val[indexCount][0];
				DenominatorAi->val[alleleCount + a1][0] += DenominatorAiBi->val[indexCount][0];
				NominatorBj->val[alleleCount + a1][0] += NominatorAjBj->val[indexCount][0];
				DenominatorBj->val[alleleCount + a1][0] += DenominatorAjBj->val[indexCount][0];
				
				NominatorAi->val[alleleCount + a2][0] += NominatorAiBi->val[indexCount][0];
				DenominatorAi->val[alleleCount + a2][0] += DenominatorAiBi->val[indexCount][0];
				NominatorBj->val[alleleCount + a2][0] += NominatorAjBj->val[indexCount][0];
				DenominatorBj->val[alleleCount + a2][0] += DenominatorAjBj->val[indexCount][0];
				
				indexCount++;
			}
		}
		alleleCount +=fl.alleles[lA];
	}
	
}


/*
	Caclulates both the gamma matrix between pairs of focal alleles Ai and Bj.
	if tolerance is >0, it computes gamma * with some tolerance
 
*/


void GammaMultiAllele_AiBj(struct landscape fl, int noAlleles, struct matrix * NominatorAiBiAjBj, struct matrix * DenominatorAiBiAjBj, struct matrix * NominatorAiBj, struct matrix * DenominatorAiBj, double tolerance)
{
	int a1,a2,lA,b1,b2,lB;
	int alleleCountA = 0, alleleCountB = 0, indexCountA = 0, indexCountB = 0;
	
	
	for(a1 = 0; a1 < noAlleles; a1++)
	{
		for(a2 = 0; a2 < noAlleles; a2++)
		{
			NominatorAiBj->val[a1][a2] = 0.;
			DenominatorAiBj->val[a1][a2] = 0.;
		}
	}
	
	for(lA = 0; lA < fl.nlocus; lA++)
	{
		for(a1 = 0; a1 < fl.alleles[lA]; a1++)
		{
			for(a2 = a1+1; a2 < fl.alleles[lA]; a2++)
			{
				indexCountB = 0;
				alleleCountB = 0;
				for(lB = 0; lB < fl.nlocus; lB++)
				{
					for(b1 = 0; b1 < fl.alleles[lB]; b1++)
					{
						for(b2 = b1+1; b2 < fl.alleles[lB]; b2++)
						{
							if(lA == lB)
							{
								NominatorAiBj->val[alleleCountA + a1][alleleCountB + b1] = 1.;
								NominatorAiBj->val[alleleCountA + a1][alleleCountB + b2] = 1.;
								
								NominatorAiBj->val[alleleCountA + a2][alleleCountB + b1] = 1.;
								NominatorAiBj->val[alleleCountA + a2][alleleCountB + b2] = 1.;
								
								DenominatorAiBj->val[alleleCountA + a1][alleleCountB + b1] = 1.;
								DenominatorAiBj->val[alleleCountA + a1][alleleCountB + b2] = 1.;
								
								DenominatorAiBj->val[alleleCountA + a2][alleleCountB + b1] = 1.;
								DenominatorAiBj->val[alleleCountA + a2][alleleCountB + b2] = 1.;
							}
							else
							{
								NominatorAiBj->val[alleleCountA + a1][alleleCountB + b1] += NominatorAiBiAjBj->val[indexCountA][indexCountB];
								NominatorAiBj->val[alleleCountA + a1][alleleCountB + b2] += NominatorAiBiAjBj->val[indexCountA][indexCountB];
								
								NominatorAiBj->val[alleleCountA + a2][alleleCountB + b1] += NominatorAiBiAjBj->val[indexCountA][indexCountB];
								NominatorAiBj->val[alleleCountA + a2][alleleCountB + b2] += NominatorAiBiAjBj->val[indexCountA][indexCountB];
								
								DenominatorAiBj->val[alleleCountA + a1][alleleCountB + b1] += DenominatorAiBiAjBj->val[indexCountA][indexCountB];
								DenominatorAiBj->val[alleleCountA + a1][alleleCountB + b2] += DenominatorAiBiAjBj->val[indexCountA][indexCountB];
								
								DenominatorAiBj->val[alleleCountA + a2][alleleCountB + b1] += DenominatorAiBiAjBj->val[indexCountA][indexCountB];
								DenominatorAiBj->val[alleleCountA + a2][alleleCountB + b2] += DenominatorAiBiAjBj->val[indexCountA][indexCountB];
								
							}
							indexCountB++;
						}
					}
				alleleCountB +=fl.alleles[lB];
				}
				indexCountA++;
			}
		}
		alleleCountA +=fl.alleles[lA];
	}
	
}

/*
	Caclulates the gamma matrix for a focal locus i and j (i.e., gamma of a mutation involving locus i on other mutations and gamma of other mutations on mutations involving locus j).
	if tolerance is >0, it computes gamma * with some tolerance
 
 */


void GammaMultiAllele_ij(struct landscape fl, int noAlleles, struct matrix * NominatorAi, struct matrix * DenominatorAi, struct matrix * NominatorBj, struct matrix * DenominatorBj, struct matrix * Nominatori, struct matrix * Denominatori, struct matrix * Nominatorj, struct matrix * Denominatorj, double tolerance)
{
	int a1,lA;
	int alleleCount = 0;
	
	for(lA = 0; lA < fl.nlocus; lA++)
	{
		Nominatori->val[lA][0] = 0.;
		Denominatori->val[lA][0] = 0.;
		Nominatorj->val[lA][0] = 0.;
		Denominatorj->val[lA][0] = 0.;
		
		for(a1 = 0; a1 < fl.alleles[lA]; a1++)
		{
			Nominatori->val[lA][0] += NominatorAi->val[alleleCount][0];
			Denominatori->val[lA][0] += DenominatorAi->val[alleleCount][0];
			Nominatorj->val[lA][0] += NominatorBj->val[alleleCount][0];
			Denominatorj->val[lA][0] += DenominatorBj->val[alleleCount][0];
			alleleCount++;
		}
	}
	
}

double GammaDistance(struct landscape fl, int distance, int tolerance)
{
	int g, gM, lA, a1, HD;
	int *diffLoci = (int *)malloc( (size_t) distance*sizeof(int) );
	int *refGeno = NULL;
	double s_j, s_ij, nominator = 0., denominator = 0.;
	
	for(g = 0; g < fl.ngenotypes; g++)
	{
		refGeno = int2genotype(fl, g, refGeno);
		for(gM = 0; gM < fl.ngenotypes; gM++)
		{
			diffLoci = genotype_diff(g, gM, &fl, &HD);
			if(HD == distance)
			{
				for(lA = 0; lA < fl.nlocus; lA++)
				{
					if( isvalueinarray(lA, diffLoci, distance) == 0 )
					{
						for(a1 = 0; a1 < fl.alleles[lA]; a1++)
						{
							if(refGeno[lA] != a1)
							{
								s_j = deltafMultiAllele(fl, g, lA, a1, tolerance);
								s_ij = deltafMultiAllele(fl, gM, lA, a1, tolerance);
								nominator += s_ij * s_j;
								denominator += pow(s_j, 2.);
							}
						}
					}
				}
				
			}
		}
	}
	
	free(diffLoci);
	return ( (nominator / denominator) );
}

/*
	Return the sum of all fitness difference squared
	Basically it is the variance
*/
double gamma_denominator(struct landscape h, double tolerance)
{

	double temp=0;
	unsigned long g;
	int j;
	
	
	for (g=0; g<h.ngenotypes; g++)
		for(j=0;j<h.nlocus;j++){
			temp += pow( deltaf(h,g,j, tolerance), 2 );
		}

	return (temp);
}


/*
	Return the sum of all fitness difference squared
	discarding values when i and j are equal
*/
double gamma_denominator_i(struct landscape h, int i, double tolerance)
{

	double temp=0;
	unsigned long g;
	int j;
	
	for(j=0;j<h.nlocus;j++)
		if(j!=i){
			for (g=0; g<h.ngenotypes; g++)
				temp+=pow(deltaf(h,g,j, tolerance),2);
		}
		
	return (temp);
}

/*
	Return the sum of all fitness difference squared
	focusing on changes at locus j
*/
double gamma_denominator_j(struct landscape h, int j, double tolerance){

	double temp=0;
	unsigned long g;


	for (g=0; g<h.ngenotypes; g++)
		temp += pow( deltaf(h,g,j, tolerance) , 2 );


	return temp;
}


/*
	The covariance of gamma for a given genotype g. See definition of eq 1
*/
double gamma_numerator(struct landscape h, int g, double tolerance){

	double temp=0;
	int l=h.nlocus;
	int i,j;
	
	for(j=0;j<l;j++) 
		for(i=0;i<l;i++)
			if(j!=i){
				temp+=deltaf(h,g,j, tolerance)*deltaf(h,gmod(h,g,i),j, tolerance);
			}
	
	if ( l == 1 )
		printf("gamma_numerator:should never happen"), exit(1);
		
	return temp/(l-1);
}
/*
	Summed over all genotypes
*/
double gamma_numerator_global(struct landscape h, double tolerance){

	double temp=0;
	int g;

	for(g=0;g<h.ngenotypes;g++)
		temp+=gamma_numerator(h,g,tolerance);

	return temp;
}



/*
	The covariance of gamma i->. See definition of eq 7
*/
double gamma_numerator_i(struct landscape h, int g, int i, double tolerance){

	double temp=0;
	int j;

	for(j=0;j<h.nlocus;j++)
		if(j!=i) {
			temp += deltaf(h,g,j, tolerance)*deltaf(h,gmod(h,g,i),j, tolerance);
		}

	return temp;
}
double gamma_numerator_i_global(struct landscape h, int i, double tolerance){

	double temp=0;
	int g;

	for(g=0;g<h.ngenotypes;g++)
		temp+=gamma_numerator_i(h,g,i,tolerance);

	return temp;
}


/*
	The covariance of gamma ->j. See definition of eq 8
*/
double gamma_numerator_j(struct landscape h, int g, int j, double tolerance){

	double temp=0;
	int i;

	for(i=0 ; i<h.nlocus ; i++)
		if(i!=j) {
			temp+=deltaf(h,g,j, tolerance)*deltaf(h,gmod(h,g,i),j, tolerance);
		}

	if ( h.nlocus == 1 )
		printf("gamma_numerator_j:should never happen"), exit(1);

	return temp/(h.nlocus-1);
}


double gamma_numerator_j_global(struct landscape h, int j, double tolerance){

	double temp=0;
	int g;
		
	for(g=0;g<h.ngenotypes;g++)
		temp+=gamma_numerator_j(h,g,j,tolerance);
	
	return temp;
}


/*
	The covariance of gamma i->j. See definition of eq 9
*/
double gamma_numerator_ij(struct landscape h, int g, int i, int j, double tolerance){

	double temp=0;
	
	if(i!=j)
		temp += deltaf(h,g,j, tolerance)*deltaf(h,gmod(h,g,i),j, tolerance);

	return temp;
}
double gamma_numerator_ij_global(struct landscape h, int i, int j, double tolerance){

	double temp=0;
	int g;
	
	for(g=0;g<h.ngenotypes;g++)
		temp+=gamma_numerator_ij(h,g,i,j,tolerance);
	
	return temp;
}


/* 
	Gamma at distance d	
*/
double gamma_numerator_dist_j_global(struct landscape h, int d, int j, double tolerance){

	double temp=0;
	int l=h.nlocus;
	int g, gf, i, v[d],
	nj[l-1];
	
	unsigned long gmax=h.ngenotypes;
	
	unsigned long indices;
	
	unsigned long binom= (unsigned long)N_choose_n(l-1,d);    /* number of neighbors at distance d, keeping locus j unchanged */


	for(i=0;i<l-1;i++)
	{
		nj[i]=(i<j)?i:i+1;   /* all locus but j */
	}
	
	
	for(i=0;i<d;i++)
		v[i]=i;   /* it contains numbers from 0 to d-1 */
	
	
	for(indices=0 ; indices<binom ; indices++){  /* for all neighbors */
	
		for(g=0;g<gmax;g++)
		{
		
			gf=g;
			
			for(i=0;i<d;i++)
			{
				gf = gmod( h, gf, nj[v[i]] );    /* generate a genome that is at distance
				                                    d and for which j is held identical changing the loci defined by v[] */
			}
			
			temp += deltaf( h,g,j, tolerance )*deltaf( h,gf,j, tolerance );
			/* check: printf("check %u : %u %u %u %f %f\n",j,g,gf,nj[v[0]]+1,deltaf(h,g,j, tolerance),deltaf(h,gf,j, tolerance)); */
		}
		
		
		/* update v[i] to explore all adequate neighbors */
		
		i=0;
		while( (v[i+1]==v[i]+1) && (i+1<d) )
		{
			i++;
		}
		
		v[i]++;
		
		for(i=i-1;i>=0;i--)
		{
			v[i]=i;
		} 
	}
	
	if (binom==0)
		printf("gamma_numerator_dist_j_global:should never happen"), exit(1);
	
	return (temp/binom);
}

double gamma_numerator_dist_global(struct landscape h, int d, double tolerance){

	double temp=0;
	int j;

	for (j=0;j<h.nlocus;j++)
	{
		temp+=gamma_numerator_dist_j_global(h,d,j,tolerance);
	}
	
	return (temp);
}


/*
	Gamma at distance 1
*/
double gamma_global(struct landscape h, double tolerance)
{

	double v=gamma_denominator(h,tolerance);
	
	if (v==0)
	{
		return (1.);
	}
	
	return (gamma_numerator_global(h,tolerance)/v);
}


double gamma_i_global(struct landscape h, int i, double tolerance){


	double v=gamma_denominator_i(h,i,tolerance);
	if (v==0)
	{
		printf("gamma_i_global:should never happen"), exit(1);
	}
	return (gamma_numerator_i_global(h,i,tolerance)/v);
}


double gamma_j_global(struct landscape h, int j, double tolerance){

	double v=gamma_denominator_j(h,j,tolerance);
	if (v==0)
	{
		printf("gamma_j_global:should never happen"), exit(1);
	}
	
	/* printf("%f / %f\n", gamma_numerator_j_global(h,j,tolerance), v); */
		
	return (gamma_numerator_j_global(h,j,tolerance)/v);
}


double gamma_ij_global(struct landscape h, int i, int j, double tolerance){

	double v=gamma_denominator_j(h,j,tolerance);
	if (v==0)
	{
		printf("gamma_ij_global:should never happen"), exit(1);
	}
	return (gamma_numerator_ij_global(h,i,j,tolerance)/v);
}


//distance greater than 1:

double gamma_dist_j_global(struct landscape h, int d, int j, double tolerance){
	double v=gamma_denominator_j(h,j,tolerance);
	if (v==0)
	{
		printf("gamma_dist_j_global:should never happen"), exit(1);
	}
	return (gamma_numerator_dist_j_global(h,d,j,tolerance)/v);
}

/*
	Gamma at distance d
*/
double gamma_dist_global(struct landscape h, int d, double tolerance){
	double v=gamma_denominator(h,tolerance);
	if (v==0)
	{
		printf("gamma_dist_global:should never happen"), exit(1);
	}
	return (gamma_numerator_dist_global(h,d,tolerance)/v);
}


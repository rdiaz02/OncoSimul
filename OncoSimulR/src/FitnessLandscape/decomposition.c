
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include  "LinearAlgebra.h"
#include  "landscape.h"
#include  "calculus.h"
#include  "vector.h"


// Calculates the number of all subsets of size maxSize
// :: NOT NEEDED ANYMORE BUT TOO NICE TO BE DELETED ::
int computeSubsetSize(struct landscape *FL, int locus, int currentSize, int maxSize)
{
	int setSize = 0;
	
	if (currentSize < maxSize)
	{
		if(locus < FL->nlocus - 1)
		{
			setSize = FL->alleles[locus] - 1;
			setSize *= computeSubsetSize(FL, locus + 1, currentSize + 1, maxSize);
			setSize += computeSubsetSize(FL, locus + 1, currentSize, maxSize);

		}
	}
	else
	{
		if (locus < FL->nlocus)
		{
			setSize = FL->alleles[locus] - 1;
			if (locus < FL->nlocus -1)
			{
				setSize += computeSubsetSize(FL, locus + 1, currentSize, maxSize);
			}
		}
	}
		
	return (setSize);
}


void createSubsets(struct landscape *FL, int locus, int allele, int currentSize, int maxSize, int * set, int **alleleArray, vector *v)
{
	int i, j;
	
	set[currentSize - 1] = alleleArray[locus][allele];
	
	if(currentSize < maxSize)
	{
		for(i=locus + 1; i < FL->nlocus; i++)
		{
			if( (i <= FL->nlocus - maxSize + currentSize) )
			{
				createSubsets(FL, i, 0, currentSize + 1, maxSize, set, alleleArray, v);
			}
		}
		if( allele < FL->alleles[locus] - 2 )
		{
			createSubsets(FL, locus, allele + 1, currentSize, maxSize, set, alleleArray, v);
		}
		 
	}
	else //if currentSize == maxSize
	{
		vector_add(v, set, maxSize);
		for (j = allele + 1; j < FL->alleles[locus] - 1; j++)
		{
			set[currentSize - 1] = alleleArray[locus][j];
			vector_add(v, set, maxSize);
		}
	}
}

void FourierTransform(struct landscape *FL, struct matrix * fourierCoefficients, double * fourierCoefficient_total, double * fourierCoefficient_epi, int * max_f)
{
	struct matrix A, inverseA, codedData, W, solutions;
	int *geno = NULL;
	int g, g2, i, j, cAlleles;
	vector subsets;
	vector_init(&subsets, FL->ngenotypes - 1);
	
	// Creates look-up table (continuously) indexing the alleles found across all loci
	
	int **alleleArray = (int **)malloc( (size_t) FL->nlocus * sizeof(int*) );
	if( !alleleArray ) fprintf(stderr, "FourierTransform: cannot allocate matrix, bye\n"), exit(1);
	
	cAlleles = 0;
	for(i = 0; i< FL->nlocus; i++)
	{
		alleleArray[i] = (int *)malloc( (size_t) FL->alleles[i] * sizeof(int) );
		if( !alleleArray[i] ) fprintf(stderr, "FourierTransform: cannot allocate matrix, bye\n"), exit(1);
		for(j = 0; j < FL->alleles[i] - 1; j++)
		{
			alleleArray[i][j] = cAlleles;
			cAlleles++;
		}
	}
	
	for(i = 1; i < FL->nlocus + 1; i++) // number of elements
	{
		for(j = 0; j < FL->nlocus - i + 1; j++)
		{
			int * set = (int *)malloc( (size_t) i * sizeof(int) );
			createSubsets(FL, j, 0, 1, i, set, alleleArray, &subsets);
			free(set);
		}
	}
	
	codedData = MemMat(FL->ngenotypes, FL->nalleles - FL->nlocus);	/* this one will contain all spin classes */
	A = MemMat(FL->ngenotypes, FL->ngenotypes);						/* this one will contain all genotypes written as spin-class products */
	inverseA = MemMat(FL->ngenotypes, FL->ngenotypes);				/* this one will contain all spin-transformed genotypes */
	W = MemMat(FL->ngenotypes, 1);									/* this one will contain the fitness */
	solutions = MemMat(FL->ngenotypes, 1);									/* this one will contain the fitness */

	
	for (g = 0; g < FL->ngenotypes; g++)
	{
		geno = int2genotype(*FL, g, geno);
		cAlleles = 0;
		for (j = 0; j < FL->nlocus; j++)
		{
			for (i = 1; i < FL->alleles[j]; i++)
			{
				if(geno[j] == i)
				{
					codedData.val[g][cAlleles] = 1.;
				}
				else
				{
					codedData.val[g][cAlleles] = -1.;
				}
				cAlleles++;
			}
		}
		W.val[g][0] = (float) FL->fitness[g];
		
		for (g2 = 0; g2 < FL->ngenotypes; g2++)
		{
			A.val[g][g2] = 1;
		}
	}
	
	for (g = 0; g < FL->ngenotypes; g++)
	{
		for (g2 = 1; g2 < FL->ngenotypes; g2++)
		{
			for(i = 0; i < vector_lengthElement(&subsets, g2-1); i++)
			{
				A.val[g][g2] *= codedData.val[g][vector_get(&subsets, g2-1)[i]];
			}
		}
	}

	
	inverseMatrix(&A, &inverseA);
	MatrixProductRef(&inverseA, &W, &solutions, 1);
	
	fourierCoefficients->val[0][0] = pow(solutions.val[0][0], 2);
	
	
	for (g = 1; g < FL->ngenotypes; g++)
	{
		fourierCoefficients->val[vector_lengthElement(&subsets, g-1)][0] += pow(solutions.val[g][0], 2);
	}
	
	
	for(i = 0; i < FL->nlocus; i++)
	{
		free(alleleArray[i]);
		
		if (i + 1 > 1)
		{
			(*fourierCoefficient_epi) += fourierCoefficients->val[i+1][0];
		}
		(*fourierCoefficient_total) += fourierCoefficients->val[i+1][0];
		(*max_f) = ( fourierCoefficients->val[i+1][0] > fourierCoefficients->val[(*max_f)][0] )?(i+1):(*max_f);
	}
	
	vector_free(&subsets);
	free(alleleArray);
	FreeMat(&A);
	FreeMat(&solutions);
	FreeMat(&codedData);

}

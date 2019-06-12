#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "calculus.h"


/*
	type could be 'i' (int), 'l' (long), 'f' (float) or 'd' (double)
*/

double *Moments(   void * untyped_array , int n, char type )
{

	double sum[3]={0,};
	double *moment;
	double val=0;
	int i,j;
	
	double epsilon=1e-6;
	
	int *iarray=NULL;
	long *larray=NULL;
	float *farray=NULL;
	double *darray=NULL;
	
	
	switch(type){
		
		case 'i':
			iarray=(int *)untyped_array;
			break;
		
		case 'l':
			larray=(long *)untyped_array;
			break;
		
		case 'f':
			farray=(float *)untyped_array;
			break;
		
		case 'd':
			darray=(double *)untyped_array;
			break;

		default:
			fprintf(stderr, "Moments: type of array incorrect, please use i, l, f or d\n");
			exit(4);
	}

	
	moment = (double *)calloc( 3, sizeof(double) );
	if(!moment)fprintf(stderr, "Moments: cannot allocate moment, bye\n"), exit(3);
	
	for(i=0;i<n;i++){
	
		switch(type){
		
			case 'i':
				val = iarray[i];
				break;
		
			case 'l':
				val = larray[i];
				break;
		
			case 'f':
				val = farray[i];
				break;
		
			case 'd':
				val = darray[i];
				break;
		}
	
		for(j=0;j<3;j++)
			sum[j] += (double )pow( val, j+1.0 );
	}
	

	/*
		mean
	*/
	moment[0] = sum[0] / n;
	
	/*
		variance
	*/
	
	if( fabs( sum[0]*sum[0] - n*sum[1] ) < epsilon )
		moment[1] = 0;
	else
		moment[1] = (sum[1]/n - pow(moment[0], 2.0));

	/*
		skewness
	*/
	if(moment[1] == 0)
		moment[2] = 0;
	else{
		moment[2] = sum[2]/n - 3.0*moment[0]*moment[1] - pow(moment[0],3);
		moment[2] /= pow( moment[1], 1.5 );
	}


	/*
		Correct for bias
	*/
	moment[1] *= ( n / (n-1.0) );
	moment[2] *= sqrt( n*(n-1.0) ) / (n-2.0); 
	

	return moment;
}

/*
	All to compute n choose n1
*/
long double LOG_factoriel_step(int n, int steps)
{

	if (steps > n)
		printf("step should be smaller than n, bye\n"), exit(2);
	if (steps == 0)
		return 0;
	return logl(n) + LOG_factoriel_step(n - 1, steps - 1);

}
long double LOG_factoriel(int n)
{

	if (n <= 1)
		return 0;
	return logl(n) + LOG_factoriel(n - 1);
}
long double log_N_choose_n(int n, int n1)
{

	long double     f = 0;
	int             n2 = n - n1;

	if (n1 > n2)
		f = LOG_factoriel_step(n, n - n1) - LOG_factoriel(n2);
	else
		f = LOG_factoriel_step(n, n - n2) - LOG_factoriel(n1);

	return f;
}

long N_choose_n(int n, int n1)
{
	return lround(exp(log_N_choose_n(n, n1)));
}






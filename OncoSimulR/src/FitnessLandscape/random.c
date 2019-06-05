/*
	random.h
	created: June 1, 2010
	
	author: gachaz
	
	function: store all functions to generate random numbers
*/


#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>


#define EPS FLT_MIN 
/*
in the original numer it was
#define EPS 1.2e-7 
*/

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define RNMX (1.0-EPS)


/***
	Random number generation
****/

static unsigned long idum_ran1;   /* the key value of ran1() */



/*
	Generate a random number between 0.00 and 1.00
	Should be called with a negative long the first time,
	then should not be reseeded.
	From Numerical recipies - ran1.c
*/
float ran1(long *idum_ran1){

	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum_ran1 <= 0 || !iy) {
		if (-(*idum_ran1) < 1) *idum_ran1=1;
		else *idum_ran1 = -(*idum_ran1);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum_ran1)/IQ;
			*idum_ran1=IA*(*idum_ran1-k*IQ)-IR*k;
			if (*idum_ran1 < 0) *idum_ran1 += IM;
			if (j < NTAB) iv[j] = *idum_ran1;
		}
		iy=iv[0];
	}
	k=(*idum_ran1)/IQ;
	*idum_ran1=IA*(*idum_ran1-k*IQ)-IR*k;
	if (*idum_ran1 < 0) *idum_ran1 += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = *idum_ran1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
/*
	Seed the function ran1() with a negative long
	Should be used only ONCE before using ran1()
*/
void seed_ran1( long seed ){

	extern unsigned long idum_ran1;            /* argument of ran1() */
	long tmp;                                  /* a tmp value fom seeding ran1();  */
	
	
	if(seed == -1)
		idum_ran1 = (unsigned long)time(NULL);     /* time_t is a long in macosx */
	else
		idum_ran1 =(unsigned long)seed;
		
	tmp = -idum_ran1;
	ran1(&tmp);                                /* seed the function */
}


/* 
	uniform_dev return a number [0,1[
*/
double uniform_dev( void ){

	extern unsigned long idum_ran1;            /* argument of ran1() */

	double dum= ran1((long *)&idum_ran1);
	
	if(dum == 0.0)
		fprintf(stderr, "BUG: ran1() can also return 0\n"),exit(1);

	if(dum == 1.0)dum=0;
	
	return dum;
	
}
/*
	uniform_int is [min , max]
*/
long uniform_int(long min, long max) 
{ 
	double dum;
 	long diff=max-min+1;
	
	dum = uniform_dev();
	 
	return (min + (int)floor(dum*(double)diff)); 
}
/*
	Get two different long in a range
*/
void uniform_2DiffInt(long min, long max, long *int1, long *int2){

	*int1 = uniform_int( min,  max);
	do{
		*int2 = uniform_int( min,  max);
	}while(*int1 == *int2);

}


/*
	Get an array of different long in a range with an exception (i.e. an excluded value);
*/
long *uniform_ArrayDiffInt(long min, long max, long n, long exception ){

	long *array,       /* returned array */
	     *array_t;;    /* a tmp array that is used to avoid redrawing the same int */

	long x,            /* a counter */
	     r;            /* random int */
	
	char e=0;
	
	
	if( exception >= min && exception <= max  )
		e=1;
	

	if(n>max-min+1-e){
		fprintf(stderr, "uniform_ArrayDiffInt: array size should be less (or equal) than the number of different values\n");
		return NULL;
	}
	
	array = (long *)malloc(  (size_t) n*sizeof(long) );
	array_t = (long *)malloc(  (size_t) (max-min+1-e)*sizeof(long) );
	if(! array || !array_t )fprintf(stderr, "uniform_ArrayDiffInt: cannot allocate array, bye\n"), exit(3);

	/*
		Fill array_t with all possible values
	*/
	e=0;
	for(x=0; x < max-min+1; x++){
	
		if( x == exception )
			e=1;
		else
			array_t[x-e] = min+x;
	}

	/*
		draw from array_t
	*/
	x=0;
	do{
		r = uniform_int( 0,  max-min-e-x );
		array[ x ] = array_t[ r ];
		array_t[ r ] = array_t[ max-min-e-x  ];
		x++;
	}while( x < n );


	free(array_t);

	return array;
}


float *uniform_array( long n  ){


	float * array, *p;
	
	array = (float *) malloc( (size_t) n*sizeof(float) );
	if(!array)fprintf(stderr,"uniform_array: cannot allocate array, bye\n"), exit(3);
	
	
	for( p=array; p-array<n; p++)
		*p=uniform_dev(  );

	return array;

}




/*
	inspired from numrec
	normal_CR_dev is mean 0 and var 1
*/
float normal_CR_dev( )
{
	static int iset=0;
	static float gset;
	float fac,rsq,v1,v2;


	if  (iset == 0) {
	
		do {
			v1=2.0*uniform_dev()-1.0;
			v2=2.0*uniform_dev()-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		
		return v2*fac;
	}
	else
	{
		iset=0;
		return gset;
	}
}

float normal_dev( double mean, double stdev ){
	return normal_CR_dev( )*stdev + mean;
}


float *normal_array( int n, double mean, double stdev ){


	float * array, *p;
	
	array = (float *) malloc( (size_t) n*sizeof(float) );
	if(!array)fprintf(stderr,"normal_array: cannot allocate array, bye\n"), exit(3);
	
	
	for( p=array; p-array<n; p++)
		*p=normal_dev( mean, stdev );

	return array;

}

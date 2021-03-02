/*
	random.h
	created: June 1, 2010
	
	author: gachaz
	
	function: header of random numbers
*/



/*
	Seed the ran1, to be done at the begining
*/
void seed_ran1_time();

/*
	Seed the ran1, to be done at the begining
*/
void seed_ran1( long seed );

/*
	Generate a random number between 0.00 and 1.00
	Should be called with a negative long the first time,
	then should not be reseeded.
	From Numerical recipies - ran1.c
*/
float ran1(long *idum_ran1);

/* 
	uniform_dev return a number [0,1[
*/
double uniform_dev( void );

/*
	uniform_int is [min , max]
*/
long uniform_int(long min, long max);

/*
	Get two different long in a range
*/
void uniform_2DiffInt(long min, long max, long *int1, long *int2);


/*
	Get an array of different long in a range with an exception (i.e. an excluded value);
*/
long *uniform_ArrayDiffInt(long min, long max, long n, long exception );





/*
	An array of random uni_dev
*/
float *uniform_array( long n  );


/*
	Normal(0,1)
*/
float normal_CR_dev( );

/*
	Normal(mu, sigma)
*/
float normal_dev( double mean, double stdev );

/*
	An array of normal_dev
*/
float *normal_array( int n, double mean, double stdev );



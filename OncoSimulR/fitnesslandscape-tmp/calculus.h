
#ifndef _MOMENTS_H_
#define _MOMENTS_H_


/*
	type could be 'i', 'l', 'f' or 'd'
*/
double *Moments(   void * untyped_array , int n, char type );

/*
	What one needs to compute N choose k
*/

long double LOG_factoriel_step(int n, int steps);
long double LOG_factoriel(int n);
long double log_N_choose_n(int n, int n1);
long N_choose_n(int n, int n1);




#endif

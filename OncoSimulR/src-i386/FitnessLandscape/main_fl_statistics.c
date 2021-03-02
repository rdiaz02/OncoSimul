
#include <stdio.h>
#include <time.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>

#include "landscape.h"
#include "genotypes.h"
#include "summary_statistics.h"

#include "calculus.h"

#include "verbose.h"

#define DATE "18 NOV 2013"



void usage( char *argv0 ){
	printf("usage is %s [options] filename\n", argv0);
	printf("version: %s\n", DATE);
	printf("to show detailed options type %s -h\n", argv0);
}


void print_help( char *argv0 ){

	usage( argv0 );

	printf("\t-R # : When comparing 2 genotypes, count them as different only if their ratio is stricly larger than #\n");
	printf("\t-z   : replace missing fitness values with 0 (otherwise check that all values are specified)\n");
	printf("\t-l   : use log fitness\n");
	printf("\t-v   : verbose --output extra stuff--\n");
	printf("\t-s   : give short output\n");
	printf("\t-o # : output is appended at the end of file #\n");
	
}



int main( int argc, char ** argv ){


	int opt;
	extern char *optarg;               /* external variables, see man 3 getopt */
	extern int optind;

	struct landscape h;                /* the ladscape itself */

	char * filename;                 /* input file */
	

	float FitnessRatio=1;             /* the ratio between two fitness should increase this value to be considered as different */
	short opt_zero=0;              /* when set to 1, replace missing fitness values by 0 */
	char opt_log=0;
	char opt_short=0;


	char *outfile=NULL;

	/* prevent confusion with the other extern verbose */
	/* char verbose_fl_stats=0; */
	verbose = 0;
	
	/*
		Parse Input
	*/

	while( (opt = getopt(argc, argv, "hR:zlsvo:")) != -1 )
	{
	
		switch(opt){
				
			case 'h':
				print_help(argv[0]);
				exit(1);
				
			case 'R':
				FitnessRatio=atof(optarg);
				break;

			case 'z':
				opt_zero=1;
				break;
				
			case 'l':
				opt_log=1;
				break;
				
			case 'v':
				/* verbose_fl_stats=1; */
				verbose=1;				
				break;
				
			case 's':
				opt_short=1;
				break;
				
			case 'o':
				outfile=optarg;
				break;
				
				
		}
		
	}

	/*
		Deal with options and syntax
	*/


	if(argc-optind != 1)usage(argv[0]),exit(1);
	filename= argv[optind+0] ;
	
	/*
		read and print landscape
	*/
	
	h = ReadFile( filename, opt_zero );
	
	/* if(verbose_fl_stats) */
	if(verbose)	  
		print_landscape( &h );
	
	OutputSummaryStats( &h, FitnessRatio, opt_log, opt_short, "coco", outfile );

	free_landscape( &h );

	return 0;
}






/*

	Generate a random FL using some rules and output it.
	Usefull to generate examples

*/
#include <stdio.h>
#include <time.h>
#include <unistd.h>
#include <stdlib.h>

#include "landscape.h"
#include "random.h"
#include "summary_statistics.h"
#include "models.h"
#include "verbose.h"



void usage( char *prog ){

	fprintf(stderr, "Version of Nov 18 2014\n");
	fprintf(stderr, "Usage is '%s [opt] {landscape} n a'\n", prog);

	fprintf(stderr, "built a fitness landscape and output it in fl format");
	fprintf(stderr, "\tOR\n");
	fprintf(stderr, "output summary statistics when option -T is given\n");

	fprintf(stderr, "Mandatory\n");
	fprintf(stderr, "\tn:  number of locus\n");
	fprintf(stderr, "\ta:  number of alleles per locus\n");
	
	fprintf(stderr, "  [landscapes]\n");
	
	fprintf(stderr, "{landscape} is defined by setting one or more of the following options\n");
	
	fprintf(stderr, "    [Fix]\n");
	fprintf(stderr, "\t-f #  : [main] add a fixed contribution\n");

	fprintf(stderr, "    [Mult]\n");
	fprintf(stderr, "\t-s #  : [main] mean fitness selection coefficient per locus ( fitness is then 1+s or s in log-scale)\n");
	fprintf(stderr, "\t-S #  : stdev for fitness for each locus (for normal_dev(); if 0 use fix fitness)\n");
	fprintf(stderr, "\t-d #  : can set a diminishing (negative) or increasing (positive) return as you approach the peak --default 0--\n");

	fprintf(stderr, "    [HoC]\n");
	fprintf(stderr, "\t-H #  : [main] stdev for fitness for House-of-Cards (normal centered on 0)\n");

	fprintf(stderr, "    [NK]\n");
	fprintf(stderr, "\t-K #  : [main] Kaufman NK with N=n and k=#\n");
	fprintf(stderr, "\t-r    : interacting loci are chosen at random --by default, they are the neighbors-- (for landscapes k or i)\n");

	fprintf(stderr, "    [Ising]\n");
	fprintf(stderr, "\t-i #  : [main] mean cost for incompatibility (if 0, use uni_dev() \n");
	fprintf(stderr, "\t-I #  : stdev cost for incompatibility (for normal_dev(); if 0 use fix cost)\n");
	fprintf(stderr, "\t-c    : last interacs with first (loci are arranged on a circle) \n");

	fprintf(stderr, "    [EggBox]\n");
	fprintf(stderr, "\t-e #  : [main] every other genotype is +/- e\n");
	fprintf(stderr, "\t-E #  : add noise on the mean effect for eggbox\n");
	
	fprintf(stderr, "    [Optimum]\n");
	fprintf(stderr, "\t-o #  : [main] optimum to target\n");
	fprintf(stderr, "\t-O #  : the spread around the optimal point, the larger the less selection for optimum\n");
	fprintf(stderr, "\t-p    : [main] the mean production value for each non 0 allele\n");
	fprintf(stderr, "\t-P    : the associated stdev\n");
	
	fprintf(stderr, "    [misc]\n");
	fprintf(stderr, "\t-a    : Specify \"a_1:..:a_i:...:a_L\", with a_i, the number of alleles at locus i \\in [1,L]\n");
	fprintf(stderr, "\t-x #  : set the seed for random generator\n");
	fprintf(stderr, "\t-n xx : set the name for FL (used in output)\n");
	fprintf(stderr, "\t-L    : All fitness are computed AND reported in log-scale\n");
	fprintf(stderr, "\t-l    : summary statistics are computed in log-scale\n");
	
	fprintf(stderr, "    [info]\n");
	fprintf(stderr, "\t-V    : be verbose\n");
	fprintf(stderr, "\t-T #  : generate # model fitness landscapes and output their statistics\n");
	
	fprintf(stderr, "\n");
}


int main( int argc, char ** argv ){


	struct landscape fl;           /* the landscape that will be implemented */
	int n,                       	    /* number of loci -- ignored when option a is set */
	    a;                       	    /* number of alleles per loci -- ignore when option a is set -- */
	
	char *allele_arg=NULL;       	    /* replace n and a when option -a is set */
	char *name="model";                 /* name of the FL --used for output */

	
	int *alleles;                       /* an array with the number of alleles at each locus */
	
	int i;

	/* extern char verbose; */                  /* level of verbosity */
	/* char verbose;     */
	
	int opt;
	extern char *optarg;                /* external variables, see man 3 getopt */
	extern int optind;
	

	struct model_opt myoptions;


	
	/* misc */
	
	short opt_log_generate=0;                /* if 1 generate FL in log scale  */
	short opt_log_stat=0;                    /* if 1 summary stat are computed after a log tranform */
	
	long init_seed=0;
	
	char opt_statistics = 0;            /* when set to 1, do r replicates of FL and output summary statistics */
	int replicates = 1;
	int r;

	verbose=0;

	myoptions = init_model( );


	/*
		Parse Input
	*/

	while( (opt = getopt(argc, argv, "s:S:H:i:I:cd:K:re:E:x:a:VT:Llo:O:p:P:f:n:")) != -1 )
	{
	
		switch(opt){

			case 'f':
				myoptions.fix = atof(optarg);
				break;
				
			case 's':
				myoptions.mu_s = atof(optarg);
				break;

			case 'S':
				myoptions.sigma_s = atof(optarg);
				break;
				
			case 'H':
				myoptions.sigma_hoc=atof(optarg);
				break;
				
			case 'i':
				myoptions.mu_ising = atof(optarg);
				break;
				
			case 'I':
				myoptions.sigma_ising = atof(optarg);
				break;

			case 'c':
				myoptions.circular_ising=1;
				break;
				

			case 'd':
				myoptions.DimRet = atof(optarg);;
				break;

				
			case 'K':
				myoptions.kaufman_K=atoi(optarg);
				break;
				
			case 'r':
				myoptions.kaufman_rand=1;
				break;
				
												
			case 'e':
				myoptions.mu_eggbox = atof(optarg);
				break;
				
			case 'E':
				myoptions.sigma_eggbox = atof(optarg);
				break;

			case 'p':
				myoptions.mu_prod=atof(optarg);
				break;
				
			case 'P':
				myoptions.sigma_prod=atof(optarg);
				break;
				
			case 'o':
				myoptions.mu_optimum=atof(optarg);
				break;
				
			case 'O':
				myoptions.sigma_optimum=atof(optarg);
				break;

			case 'x':
				init_seed=atol(optarg);
				break;
				
			
			case 'a':
				allele_arg=optarg;
				break;
				

			case 'V':
				verbose=1;
				break;
				
			case 'T':
				opt_statistics=1;
				replicates = atol(optarg);
				break;
				
			case 'L':
				opt_log_generate=1;				
				break;
				
			case 'l':
				opt_log_stat=1;				
				break;
				
			case 'n':
				name=optarg;				
				break;
				
		}
		
	}


	/*
		Deal with options and syntax
	*/

	if(argc-optind != 2 )usage(argv[0]),exit(1);

	n=atoi( argv[optind+0] );
	a=atoi( argv[optind+1] );

	/*
		Retrieve alleles or built the array
	*/
	if( allele_arg )
		alleles=char_arg2array_int(  allele_arg , &n );
	else
	{
		alleles=(int *)malloc( (size_t) n*sizeof(int) );
		if(!alleles)fprintf(stderr, "cannot allocate memory for alleles, bye\n"), exit(1);
		for(i=0;i<n;i++)
			alleles[i]=a;
	}
	//free(alleles);
	
	/*
		See ran
	*/
	if(init_seed == 0)
		seed_ran1( (long) time(NULL) );
	else
		seed_ran1( init_seed );
	
	
	/*
		Init
	*/
	init_landscape( &fl, n, alleles );

	
	/*
		Built FL
	*/
	for( r=0; r<replicates ; r++ ){
	
	
		MixedModel( &fl, myoptions, opt_log_generate );         /* if opt_log_generate == 1, all fitness values are computed and reported in log-scale */
		
		if(opt_statistics == 0)
			output_landscape( &fl );
		else
		{
			OutputSummaryStats( &fl, 1, opt_log_stat , (r==0)?1:2, name, NULL ); 
		}

	}

	free(alleles);
	free_landscape( &fl );
	
	
	return 0;
}






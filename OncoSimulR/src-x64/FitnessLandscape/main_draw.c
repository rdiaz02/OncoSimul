/*

        Draw FL from data. Usefull to get a visual sense of the landscape

	

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>


#include "landscape.h"
#include "genotypes.h"
#include "drawings.h"



/*
	Help
*/
void printf_help()
{
	printf("Display artistic landscapes. options are:\n");

	printf("\t-h   : this 'h'elp\n");

	printf("\n  [global]\n");
	printf("\t-S # : Change the Scale (height and width) by a factor (default is 1.0)\n");
	printf("\t-H # : Change the Height by a factor (default is 1.0)\n");
	printf("\t-R # : place that 'R'eference on left\n");
	printf("\t-s   : draw compact genotype\n");
	printf("\t-f   : draw flat landscape\n");

	printf("\n  [fitness]\n");
	printf("\t-z   : replace missing fitness values with 0 (otherwise treat them as unknown)\n");
	printf("\t-l   : use 'l'og of fitness values given in the input F.L.\n");
	
	printf("\n  [edges]\n");
	printf("\t-t # : 't'hreshold value for the fitness ratio to replace a plain with a dashed line (default is 2.0)\n");
	printf("\t-c   : 'c'lear out edges with fitness ratio under the 't'hreshold (do not represent dashed lines) \n");
	printf("\t-C   : draw only 'C'hains \n");
	printf("\t-r # : draw edges from a starting 'r'eference top-hill\n");

}


/*
	Automatic draws of input F.L.
*/
int main(int argc,char **argv)
{

	struct landscape land;       /* the fitness landscape itself */
	
	char *file_out,              /* input and output files */
	     *file_in;
	     
	float threshold=1.0;         /* change aspect of connection if ratio of fitness is greater than threshold */
	
	int opt_log=0,               /* consider the log(fitness) instead of the fitness */ 
	    opt_clear=0,             /* if 1, do not repsents connection with ratio of fitness < threshold */
	    opt_zero=0,              /* if 1, missing vqlues qre filled with fitness 0, otherwise, they are set to DEFAULT_FITNESS */
	    opt_ref=-1,              /* if >= 0, use it as a reference to draw connections */
	    opt_to=-1,              /* if >= 0, use it as a reference to draw connections to*/	
	    opt_compact=0, 			/* if option -s then compact genotypes instead of big*/
	    opt_flat=0,				/* if 1 flat landscpae*/	
	    opt_chains=0;           /* chains */	

	float rescale_height=1.0,
	      rescale_width=1.0;

	int GenoRef=0;               /* the genotype that is display on left */

	int R;                /* the radius -- it scales the whole picture */
	    
	int opt;                     /* to parse the options */
	extern int optind;
	extern char *optarg;

	int Only_Mutation=-1;


	while( (opt = getopt(argc, argv, "hlcfszt:r:R:H:S:M:C")) != -1 )
	{
		switch(opt){
		
			case 'h':                           /* the fabulous help */
				printf_help();
				exit(1);
		
			case 't':
				threshold = atof(optarg);    /* value of threshold */
				break;
		
			case 'r':
				opt_ref = atoi(optarg);         /* value of ref to draw connections */
				break;
		
			case 'R':
				GenoRef = atoi(optarg);         /* set the genotype that is on the left */
				break;
		
			case 'l': 
				opt_log=1;                       /* when set, use log of input values for scale */
				break;
		
			case 'c':
				opt_clear=1;	                /* when set fitness difference below threshold are not drawn */
				break;
			
			case 'C':
				opt_chains=1;	                /* just draw chains instead of full landscapes */
				break;
			
			case 'z':
				opt_zero=1;	                /* when set to 1, fill missing values with 0 */
				break;
			
			case 'H':
				rescale_height=atof(optarg);	/* used to rescale the picture on height */
				break;
				
			case 'S':
				rescale_width=atof(optarg);	/* used to rescale the picture on width */
				break;
			case 's':
				opt_compact=1;
				break;
				
			case 'f':
				opt_flat=1;
				break;


			case 'M':
				Only_Mutation=atoi(optarg);	/* used to rescale the picture on width */
				break;
		}
	}
	
	
	argc -= optind;           /* decrease argc from the number of option you get */
	if( argc != 1 )
	
		fprintf(stderr, "usage is: '%s [opt] <file>'\ntype '%s -h' for more help\n",argv[0],argv[0]),exit(1);

	
	argv += optind;         /* increase the pointer with optind */	
	file_in = *argv;


	R = 7*rescale_width;


	/*
		Read FLL and fill landscape structure
	*/
	land= ReadFile(file_in, opt_zero);
	

	/*
		Built outfile name
	*/
	file_out=malloc(sizeof(char)*strlen(file_in)+5);
	if(!file_out)fprintf(stderr, "main: cannot allocate file_out, bye"), exit(3);
	sprintf(file_out,"%s.svg",file_in);              


	/* 
		Draw the representation 1sr arg is 1 for web only
		
	*///last 0 to be changed when parameter flat is on
	draw_FL( 0, &land, file_out, threshold, opt_log, opt_clear, opt_ref,opt_to, GenoRef, rescale_height,R, Only_Mutation,0,opt_chains,1000,1000,(int)opt_compact,opt_flat,"");


	fprintf(stderr, "SVG file is : %s.svg - it can be open from a web browser (or Inkscape)\n", file_in);
       
	free(file_out);

	return 0;
}


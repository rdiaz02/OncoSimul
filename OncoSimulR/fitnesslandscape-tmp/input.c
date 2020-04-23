#include <stdio.h>
#include <stdlib.h>
#include "landscape.h"
#include <string.h>


/*
	char_arg should be numbers separated by ':'
	return an array with those numbers
*/
int * char_arg2array_int(  char *char_arg , int *nlocus ){

	
	char *p;
	int *ni;
	int i=0;
	
	*nlocus=1;
	p = char_arg;
	
	while( *(p++) )
		if( *p == ':' )
			(*nlocus)++;

	ni = (int *)calloc( (size_t)(2*(*nlocus)-1), (size_t)sizeof(int) );
	if( !ni )fprintf(stderr, "parse_sparg: cannot allocate ni, bye\n"),exit(1);

	for (p = strtok(char_arg, ":"); p; p= strtok(NULL, ":"))
		ni[i++] = atoi(p);
	
	p=char_arg;
	for(i=0;i<(*nlocus)-1;i++){
		while(*(++p));
		*p=':';
	}
	
	
	return ni;
}





int * GetOrderFromFile( char *filename ){

	FILE *f;
	char letter;
	int nlocus=1,
	    l,         /* locus */
	    g,         /* genotype */
		ch;
	int *alleles;
		
	int line=0;
	int *order;
	float fitness;
	struct landscape fl;

	f=fopen( filename, "r"  );
	if( !f )fprintf(stderr, "error while opening file %s. Make sure it exists. Bye\n", filename), exit(2);
	
	
	while( (letter = (char)fgetc(f) ) != '\n' ){
		if( letter == ' ' )
			nlocus++;
	}
	rewind(f);
	
	/*fprintf(stderr, "from file %s, there are %d locus\n", filename, nlocus);*/

	alleles = (int *) malloc( (size_t) nlocus * sizeof( int )  );
	if(! alleles )fprintf(stderr, "ReadFile: cannot allocate alleles, bye\n"), exit(3);


	for( l=0; l<nlocus; l++)
	  // fscanf( f, "%d", alleles+l );
	  if(fscanf( f, "%d", alleles+l )){};
	// the ignoring return value of ‘fscanf’ warning	

	init_landscape( &fl, nlocus, alleles);
	
	order = (int *)malloc(fl.ngenotypes*sizeof(int));
	if(!order)fprintf(stderr, "cannor allocate order, bye\n"),exit(3);
	
	line = 1;

	while( fgetc(f) != EOF )
	{
	
		fseek(f, -1, SEEK_CUR);
		
		for( l=0; l<nlocus; l++)
		{
		  // fscanf( f, "%d", alleles+l );
		  if(fscanf( f, "%d", alleles+l )){};
		}
		// fscanf( f, "%f", &fitness );
		if(fscanf( f, "%f", &fitness )){};
		
		while( ((ch = fgetc(f)) != '\n') && (ch != EOF) );
		
		line ++;

		g = genotype2int( fl, alleles );
		
		order[g]= line;

	}
	
	free_landscape( &fl );
	
	fclose(f);

	return order;
}





struct landscape ReadFile( char *filename, int opt_zero ){

	FILE *f;
	char letter;
	int nlocus=1,
	    l,        /* locus */
	    g,         /* genotype */
		ch;			/* fgetc return */

	struct landscape fl;
	
	long pos;
	
	int *alleles;
	
	float fitness;
	
	int line=0;
	

	f=fopen( filename, "r"  );
	if( !f )fprintf(stderr, "error while opening file %s. Make sure it exists. Bye\n", filename), exit(2);
	pos=ftell(f);
	letter=(char)fgetc(f);
	if (letter=='#') // first line can be a comment describing landscape
		{
	
		while( (letter = (char)fgetc(f) ) != '\n' );
		pos=ftell(f);
		}
	else
		ungetc(letter,f);
	while( (letter = (char)fgetc(f) ) != '\n' ){
		if( letter == ' ' )
			nlocus++;
	}
	//rewind(f);
	fseek(f,pos, SEEK_SET);
	/*fprintf(stderr, "from file %s, there are %d locus\n", filename, nlocus);*/

	alleles = (int *) malloc( (size_t) nlocus * sizeof( int )  );
	if(! alleles )fprintf(stderr, "ReadFile: cannot allocate alleles, bye\n"), exit(3);

	
	for( l=0; l<nlocus; l++)
	{
	  // fscanf( f, "%d", alleles+l );
		if(fscanf( f, "%d", alleles+l )){};
	}
	
	/*
	for( l=0; l<nlocus; l++)
		fprintf( stderr, "%d ", alleles[l] );
	fprintf(stderr, "\n");
	*/
	
	init_landscape( &fl, nlocus, alleles);
	
	for (g=0;g<fl.ngenotypes;g++)
		if( opt_zero )
			fl.fitness[g] = 0;
		else
			fl.fitness[g] = DEFAULT_FITNESS;
		
	line = 1;

	while( fgetc(f) != EOF )
	{
	
		fseek(f, -1, SEEK_CUR);
		
		for( l=0; l<nlocus; l++)
		{
		  //fscanf( f, "%d", alleles+l );
			if(fscanf( f, "%d", alleles+l )){};
		}
		//fscanf( f, "%f", &fitness );
		if(fscanf( f, "%f", &fitness )){};
		
		
		while( ((ch = fgetc(f)) != '\n') && (ch != EOF) );
		
		
		line ++;

		for( l=0; l<nlocus; l++ )
		{
			if( alleles[l] > fl.alleles[l]-1  )
			{
				fprintf(stderr, "line %d, allele %d has a type %d should be <= %d, bye\n", line, l, alleles[l], fl.alleles[l]-1 ), exit(4);
			}
		}
		g = genotype2int( fl, alleles );
		fl.fitness[g]= fitness;

	}

	if( line-1 != fl.ngenotypes )
		fprintf(stderr, "Warning: there are %d missing genotypes (theoretical total: %d)\n", line-1, fl.ngenotypes );
	
	/*
		Set maxf et minf --very usefull for drawing--
	*/
	g=0;
	while( fl.fitness[g] == DEFAULT_FITNESS )
		g++;
	
	fl.minf = fl.maxf = fl.fitness[g];
	
	for ( ; g<fl.ngenotypes ; g++){
	
		if( fl.fitness[g] != DEFAULT_FITNESS && fl.fitness[g]>fl.maxf )
			fl.maxf = fl.fitness[g];
		
		if( fl.fitness[g] != DEFAULT_FITNESS && fl.fitness[g]<fl.minf )
			fl.minf = fl.fitness[g];
	}
	
	/*
		Set the number of neighbors (only valid for all genotypes in a full F.L.)
	*/
	fl.neighbors=0;
	for( l=0; l<fl.nlocus; l++)
		fl.neighbors += fl.alleles[l] - 1;
	
	fclose(f);

	return fl;
}



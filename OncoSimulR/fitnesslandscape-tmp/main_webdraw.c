#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <dirent.h>
#include <sys/stat.h>
#include <time.h>
#include <float.h>
#include "landscape.h"
#include "genotypes.h"
#include "drawings.h"
#include "random.h"

#include "summary_statistics.h"
#include "models.h"
#include "web/webparam.h"
#include "web/web_common.h"
#include "web/server.h"




void html_error(int c)
{
	switch(c){
		default: 
		printf("an error %d was detected\n",c);exit(1);
		}
}
/*------------------------------------------------------*/
/*search for delimiter beetween different dorm values*/
int search_delim(char *cgiinput, char *delim)
	{
	int i=0,lequel=0;
	char *begin,*end;
	char *achercher[2]={"WebKitFormBoundary","-----------------------------"};// delimitors are parts of either first either 2nd string 

	for (i=0;i<2;i++)
		{
		begin=strstr(cgiinput,achercher[i]); // find begining of delimitor
		if (begin!=NULL)
			break;
		}
	lequel=i;	
	if (i==2)
		{
		printf("<BR>%s<BR>",cgiinput);
		html_error( 4 );
		}
 
	end=begin+strlen(achercher[lequel]); //look for delimitor 's end
	i=0;
	
	while ((*(end+i) !=' ') && (*(end+i)!=13)&& (*(end+i)!=10) && (*(end+i) !='\0')) 
			i++;
	if (*(end+i) =='\0' ||i>=256) {printf("delimiter not found %d %.128s <BR>\n",i,end);html_error( 5);} //too much char in the delim
	
	// copy the clean delimitor
	strncpy(delim,begin,strlen(achercher[lequel])+i);
	delim[strlen(achercher[lequel])+i]='\0';
	// count how many fields
	for (begin=cgiinput,end=begin,i=0;end!=NULL;i++,end=strstr(begin,delim),begin=end+1)
		; //do nothing but put the ; here because of the warning
	i--;      
	return(i-1);//last delim is for ending the form
	}
/*------------------------------------------------------*/
void printweb_landscape( struct landscape *l ){

	int g;
	int *geno=NULL;

	printf("landscape (%d loci; %d genotypes): <BR>", l->nlocus, l->ngenotypes);
	print_genotype( l->alleles, l->nlocus );
	printf("\n");

	for(g=0; g<l->ngenotypes; g++ ){
	
		geno = int2genotype( *l, g , geno );
		printf("genotype: %d '", g);
		print_genotype( geno, l->nlocus );
		printf("'  %f<BR>\n",  l->fitness[g] );

	}

	free(geno);
}

/*------------------------------------------------------*/
/** Read the CGI input and place all name/val pairs into list.        **/
/** Returns list containing name1, value1, name2, value2, ... , NULL  **/
char **getMultPartData(int *nb,char *nf) {
   	int i,j,k;
    int taille_chps; //,taille_delim;
    int content_length;
    char *cgiinput ;
    char **cgivars ;
    int paircount ;
	char delim[256];
	char *newName;
	char *begin,*end;
	

    

    // perform some validity tests on cgi input 
    if ( !(content_length = atoi(getenv("CONTENT_LENGTH"))) ) html_error( 1 );
    if ( !(cgiinput= (char *) malloc(content_length+1)) )html_error( 2 ); 
    if (!fread(cgiinput, content_length, 1, stdin)) html_error( 3 );

	cgiinput[content_length]='\0' ; 

//	printf("%s\n",cgiinput); exit(1);
	
	paircount=search_delim(cgiinput, delim);/*search the delimiter and how many of them*/
    cgivars= (char * *) malloc(2*paircount*sizeof(char *)) ;
//    taille_delim=strlen(delim);

	begin=strstr(cgiinput,delim);

	/*this part read and splits the string into an array cgivars organized as followed:
	  cgivars[i] contains the name of the field cgiVars[i+1] contains the value(s)*/
	for(i=0;i<(2*paircount);i=i+2)
		{

		end=strstr(begin+1,delim); //looking for next delim

		if(end==NULL || strlen(end)-2 <=strlen(delim) )
			break;
		taille_chps=end-begin;
		newName=strchr(begin,'=')+1;// looking for begining of name
		
		if (*newName=='"')
			newName++;
			
		
		k=0;
		while (*(newName+k) !='"' && (newName+k)!=NULL) k++; // look for end of name
		
		if ((newName+k)==NULL)html_error( 7 );
		
		if (k==0) printf("arrrg %d\n",i); //should never arrives exept if someone put a field with no name in the html form...
		cgivars[i]=malloc (sizeof(char)*(k+1));
		strncpy(cgivars[i],newName,k);
		
		cgivars[i][k]='\0';  
		

		j=0;

		cgivars[i+1]=malloc (taille_chps- k +10);
		newName=newName+k;
		if (*newName=='"')
		newName++;

		while (newName!=end && newName!=0)
			{
			
			cgivars[i+1][j++]=*(newName);
			newName++;
			if (j>taille_chps- k+6) break;
			}
		if ((j>=2) && ((cgivars[i+1][j-1]=='-' &&	cgivars[i+1][j-2]=='-'))) //on enleve les -
			{
			while(cgivars[i+1][j-1]=='-'){cgivars[i+1][j-1]='\0';j--;}
			}	
		cgivars[i+1][j]='\0';
//		printf("%d %s<BR>",i+1,cgivars[i+1]);
		begin=end +1;	
//		printf("<HR>\n");

	}	

	free(cgiinput) ;

	*nb=i/2;

    return cgivars ;
    
}

/*get the value of the field 'to_find'*/
float get_value(char *to_find,char **params,int nb_couples)
{
	int i;
	for (i=0;i<nb_couples*2;i=i+2)
		if (strcmp(to_find,params[i])==0)
			return(atof(params[i+1]));
	return(-1);		
}
int get_intvalue(char *to_find,char **params,int nb_couples)
{
	int i;
	for (i=0;i<nb_couples*2;i=i+2)
		if (strcmp(to_find,params[i])==0)
			return(atoi(params[i+1]));
	return(-1);		
}


/*get the string associated to the field 'to_find'*/
char *cleanrc(char *str)
{
	char *l;
	int i,t=strlen(str),n=0;
	
	for (i=0;i<t;i++)
		if (str[i]!=10 && str[i]!=13)
			n++;
	
	l=malloc(sizeof(char)*(n+1));	

	n=0;
	for (i=0;i<t;i++)
		if (str[i]!=10 && str[i]!=13)
			l[n++]=str[i];
	
	l[n]='\0';
	
	return l;
}


char *get_string(char *to_find,char **params,int nb_couples)
{
	int i;
	for (i=0;i<nb_couples*2;i=i+2)
		if (strcmp(to_find,params[i])==0)
			return(cleanrc(params[i+1]));
	return(NULL);		
}


/*look if 'to_find' is in the form _ used for checked box_*/
int presence(char *to_find,char **params,int nb_couples)
{
int i;
for (i=0;i<nb_couples*2;i=i+2)
	if (strcmp(to_find,params[i])==0)
		return(1);
return(0);		
}
		       
/*void draw_end()
{
printf(  "</html>\n");
printf(  "<!-- end -->\n");
}
*/


int compare_geno(int *mask,int *geno,  int nbal)
{
int i;
for (i=0;i<nbal;i++)
	if (mask[i]!=-1 && mask[i]!=geno[i])
		return 0;
return 1; //if compatible		
}



struct landscape ReadString( char *data, int opt_zero )
{
	char letter;
	int nlocus=1,
		length,
	    l,        /* locus */
	    g;         /* genotype */
	char *val = NULL;
	struct landscape fl;
	char leslocus[132];
	int *alleles;
	
	float fitness;
	
	int line=0,i=0,j=0;
	

	length	=strlen(data);
	val=malloc(sizeof(char)*(length+1));//max for val
//	 printf( "<BR>");	
		
	while( *data == 10 || *data == 13  || *data == ' '|| *data=='\t') data++ ; //1st char must be a value
	if (*data=='#') //first line is a comment 
		while( *data !=10 && *data != 13 ) data++ ;//do it again
	
	if (*data ==10 || *data == 13) data++;	

	j=i=0;
	
	//pb des espaces en plus ou des tab....
	while( (letter = (char)data[i++]) != 10 && letter!=13){ 
		
		leslocus[j++]=letter;
//		printf("%c ",letter);
		if (j>132)	
		{printf("problem in file or too much allels<BR>");exit(1);}
//		if( letter == ' ' || letter=='\t')
//			nlocus++;
//			printf("%d <br>", data[i-1]);
	}
	
	leslocus[j]='\0';	


	nlocus=trim (leslocus); //enleve espaces et renvoie le nbre de locus
	
	

 	
	alleles = (int *) malloc( (size_t) nlocus * sizeof( int )  );
	if(! alleles )fprintf(stderr, "ReadFile: cannot allocate alleles, bye\n"), exit(3);

	
	for( l=0; l<nlocus; l++)
		{
		sscanf(data, "%s", val );
		*(alleles+l)=atoi(val);

		data+=strlen(val);
		while (((*data==' ')||(*data=='\t')) && *data !='\n') data++;
		
	}

//printf("alleles read<BR>");		

	init_landscape( &fl, nlocus, alleles);
	
	for (g=0;g<fl.ngenotypes;g++)
		if( opt_zero )
			fl.fitness[g] = 0;
		else
			fl.fitness[g] = DEFAULT_FITNESS;
		
	line = 1;

while( *data == 10 || *data == 13  || *data == ' '|| *data=='\t') data++ ; //1st char must be a value
	
	
	while(*data!='\0'){
		
	
		for( l=0; l<nlocus; l++)
			{
			sscanf(data, "%s", val );
		
			*(alleles+l)=atoi(val);
			data+=strlen(val);
			if (strlen (val)>length){printf("problem in file or no rc found<BR>");exit(1);}
			while ((*data==' ' || *data=='\t') && *data !='\n') data++;
		
			}
		
			sscanf(data, "%s", val );
			fitness=strtod(val, (char **)NULL);
				
			if (strlen (val)>1320){printf("problem in file or no rc found<BR>");exit(1);}
				
			data+=strlen(val);
		
				
while( *data == 10 || *data == 13  || *data == ' '|| *data=='\t') data++ ; //1st char must be a value
	
		if (*data !='\0')
		line ++;
	
		for( l=0; l<nlocus; l++ )
			if( alleles[l] > fl.alleles[l]-1  )
				
				printf("line %d, allele %d has a type %d should be <= %d, bye\n", line, l, alleles[l], fl.alleles[l]-1 ), exit(4);

		g = genotype2int( fl, alleles );
		fl.fitness[g]= fitness;
	
	}

	
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
	
		free(val);

	return fl;
}










char *skipHeader(char *data, char **filen){	
	
	int i=0;

	char *filename,*bout;
	
	
	/*first read the filename*/
	if ((bout=strstr(data,"filename="))==NULL)
		printf("bad filename keyword was not found in <BR>%s<BR>\n",data),exit(1);
	else
		{
		i=0;
		
		data=bout+strlen("filename=")+1;
		if (*data=='"') data++;	
				
		while (*(data+i++) !='"' && *data!='\0');
		
		if (*data=='\0' ||i==0)
			printf("bad filename <BR>\n"),exit(1);
			
		filename=malloc(sizeof(char)*(i+1));
		if (filename==NULL) printf("pbmemorry<BR>"),exit(1);
		
		i=0;
		while (*(data+i)!='"')
			  {filename[i]=*(data+i);i++;}
		filename[i]='\0';	  	
//		printf("<BR>filename: %s<BR>",filename);
		}

	/* then  find where the data begins and *so funny* this differs from one navigator to another :~/ */
	bout=strstr(data,"octet-stream");
	if(bout==NULL)
		{
		bout=strstr(data,"x-fluid");
		if (bout==NULL)
			{
			bout=strstr(data,"text/plain");
			if (bout!=NULL)
				 data=bout+strlen("text/plain");
			else
				printf("Bad data or mime type (try to change the extension of your file by .txt)<BR>\n"),exit(1);
			}
		else
			data=bout+strlen("x-fluid");
		}
	else	
		data=bout+strlen("octet-stream");
	
	/*clean the end of string if some separator's chars remain*/
	i=strlen(data)-1;
	while (data[i]=='-' || data[i]==10 || data[i]==13 || data[i]==' ') data[i--]='\0';
	
	*filen=filename;

	return(data);

	}

	
int getmask(char **cgivars,struct landscape land,int *mm,int nb)
{
int i,masked=0;
char suffix[32];

for (i=0;i<land.nlocus;i++)
	{
	sprintf(suffix,"LAND_MASK_%d",i);
	mm[i]=get_value(suffix,cgivars,nb);
	
	if (mm[i]!=-1)
		masked=1;
	}

return(masked);
}



char *clean (char *data)	
{
	int i=strlen(data)-1;
	while (data[i]=='-' ||data[i]==' '||data[i]==10 ||data[i]==13) data[i--]='\0';
	i=0;
	while (data[i]=='-' ||data[i]==' '||data[i]==10 ||data[i]==13) i++;;
	return(data+i);
}


/*******************************
clean regularly temp dir
*******************************/
void menage(char *ledir)
{
	struct dirent *dp;
	DIR *md;
	time_t tps_local;
	time_t tps_fich;
	struct stat buf;
	char commande[256];
	char lepetitdir[256];
	
	md=opendir( ledir );
	tps_local=time(NULL);
	
	 while ((dp=readdir(md))!=NULL)
		   {
			
			sprintf(lepetitdir,"%s%s",ledir,dp->d_name);
			stat(lepetitdir,&buf);
			tps_fich= buf.st_mtime;
//			printf("%d fichier=%d %d %d\n",tps_local,tps_fich,tps_local-tps_fich,(float)(tps_local-tps_fich)/86400.0);		
			if ((tps_local-tps_fich) > 86400*3 && strcmp(dp->d_name,".")!= 0 && strcmp(dp->d_name,"..")!= 0)
				{ 
				
				sprintf(commande, "rm -rf %s",lepetitdir);
//				printf("%s\n",commande);			
				system(commande);
				}
//			else 
//			printf("%s %d \n",lepetitdir,difftime(tps_local,tps_fich));		
			 }
			 
		  
	 closedir(md);
}
/*
void setminMax(struct landscape *l)
{
int g;
float min=FLT_MAX;
float max=FLT_MIN;
for(g=0; g<l->ngenotypes; g++ )
	{
	if (l->fitness[g]>max)
		max=l->fitness[g];
	if (l->fitness[g]<min)
		min=l->fitness[g];
		
	}
l->minf=min;
l->maxf=max;	
}

*/

//check clean string
char * cleanRC(char *s1)
{
int i,l=strlen(s1),j;

char *n;
for (i=0,j=0;i<l;i++)
	if (s1[i]!=10 &&s1[i]!=13)
		j++;
		
n=malloc(sizeof(char ) *(j+1));

for (i=0,j=0;i<l;i++)
	if (s1[i]!=10 && s1[i]!=13)
		n[j++]=s1[i];
return n;
		
}



int main (int argc, char **argv)
{
	char **cgivars,*data, *thestats=NULL;

	char dataFilename[1000];
	int nbFields;
	int ii,ifile=-1,iprevious=-1;
	struct landscape land;
	struct landscape newland;
	struct landscape *p_land;
	int opt_masked=0;
	char *stat_formated;    
	int i,flag=0,incomp=0;
	float R;
	char *filename;
	char *file_out=NULL;
	char *file_csv=NULL;
		char *file_txt=NULL;
	struct model_opt junk;
	AllParams param;
	char postfix[1024];
	

	seed_ran1( (long) time(NULL) );	
	printf("Content-type: text/HTML\n\n");
	initparam(&param);

	cgivars=getMultPartData(&nbFields,dataFilename);
	
/*leave this commented to check easily the parsing of web page*/
//		for (ii=0;ii<nbFields*2;ii++)
//		printf("%s <BR>\n",cgivars[ii]);


	for (ii=0;ii<nbFields*2;ii+=2)
			{
			if (strcmp(cgivars[ii],"LAND_FILE")==0)
				ifile=ii;
			else
			if 	(strcmp(cgivars[ii],"LAND_THRESHOLD")==0)
	//				threshold=atof(cgivars[ii+1]);
				param.LAND_THRESHOLD=atof(cgivars[ii+1]);
						
			else	
			if 	(strcmp(cgivars[ii],"LAND_REFERENCE")==0)
	//					GenoRef=atoi(cgivars[ii+1]);
				param.LAND_REFERENCE=atoi(cgivars[ii+1]);
							
			else	
			if 	(strcmp(cgivars[ii],"LAND_DRAWFROM")==0)
	//						opt_ref=atoi(cgivars[ii+1]);
				param.LAND_DRAWFROM=atoi(cgivars[ii+1]);	
			else	
					
			if 	(strcmp(cgivars[ii],"LAND_OPTCLEAN")==0)
		//						opt_clear=1;
				param.LAND_OPTCLEAN=1;		
			else
			if 	(strcmp(cgivars[ii],"LAND_LOG")==0)
		//							opt_log=1;
				param.LAND_LOG=1;
			else	
							
			if 	(strcmp(cgivars[ii],"LAND_MISSING")==0)
		//								opt_zero=1;
		 		param.LAND_MISSING=1;
			else			
			if (strcmp(cgivars[ii],"LAND_SCALE_W")==0)
		//									rescale_width=atof(cgivars[ii+1]);
				{
				param.LAND_SCALE_W=atof(cgivars[ii+1]);
				if (param.LAND_SCALE_W==0)
				param.LAND_SCALE_W=1;
				}
			else	
			if (strcmp(cgivars[ii],"LAND_SCALE_H")==0)
		//										rescale_height=atof(cgivars[ii+1]);
				param.LAND_SCALE_H=atof(cgivars[ii+1]);
			else
			if (strcmp(cgivars[ii],"LAND_PREVIOUS_FILE")==0)
				iprevious=ii;
//				param.LAND_PREVIOUS_FILE=ii;
			else
			if 	(strcmp(cgivars[ii],"LAND_DRAWTOEND")==0)
		//											opt_to=atoi(cgivars[ii+1]);	
				param.LAND_DRAWTOEND=atof(cgivars[ii+1]);
			else
			if (strcmp(cgivars[ii],"LAND_CHAINS")==0)
		//												opt_chains=1;	
				param.LAND_CHAINS=1;
			else	
			if (strcmp(cgivars[ii],"LAND_WINDOW_WIDTHSIZE")==0)
				param.LAND_WINDOW_WIDTHSIZE=atoi(cgivars[ii+1]);
			else	
			if (strcmp(cgivars[ii],"LAND_WINDOW_HEIGHTSIZE")==0)
				param.LAND_WINDOW_HEIGHTSIZE=atoi(cgivars[ii+1]);
			else	
			if (strcmp(cgivars[ii],"LAND_ONLY_MUT")==0)
				param.LAND_ONLY_MUT=atoi(cgivars[ii+1]);
			else
			if (strcmp(cgivars[ii],"LAND_COMPACT")==0)
				param.LAND_COMPACT=1;			
			else
			if (strcmp(cgivars[ii],"LAND_FLAT")==0)
				param.LAND_FLAT=1;	
			else													
				if (strcmp(cgivars[ii],"LAND_CHECK_SAMENBR")==0) //if found then checked so sameNumber
					param.LAND_CHECK_SAMENBR=1;							
			}	
				

		
		R = 7*(param.LAND_SCALE_W/100.0);
		
		
  		if (iprevious== -1)
  		{
//  		printf("pas de previous-------->%s<BR>",cgivars[ifile+1]);
		data=skipHeader(cgivars[ifile+1], &filename);
		}
		else
		{
        filename=clean (cgivars[ifile+1]);
         data=cgivars[iprevious+1];
  //   printf("-------->%s<BR>",data);
		}
        param.LAND_MODEL=malloc(sizeof(char)*(strlen(filename)+1));    //copy file name instaed of model name    
		strcpy(param.LAND_MODEL,filename);
		
		land = ReadString(data,param.LAND_MISSING);
		
		param.LAND_MASK=malloc(sizeof(int)*land.nlocus);
		
		for (i=0;i<land.nlocus;i++)
			param.LAND_MASK[i]=-1;
		opt_masked= getmask(cgivars,land,param.LAND_MASK,nbFields);	 
	if (land.minf< 0 && param.LAND_LOG)
		{printf("<BR>Cant use log when negative fitness value. Please unselect log and re run<BR>");exit(1);}
		sprintf(postfix,"%d.%ld",getpid(),random());

	file_out=malloc(sizeof(char)* (strlen(filename)+6+strlen(OUTDIR)+strlen(postfix))); //6= 2points  and svg will be added in filename below
	file_csv=malloc(sizeof(char)* (strlen(filename)+6+strlen(OUTDIR)+strlen(postfix)));
	if(!file_out)fprintf(stderr, "main: cannot allocate file_out, bye"), exit(3);
		srand(time(NULL));
	//create dir where files are going to be written
	
	sprintf(file_out,"%s%s.%s.svg",OUTDIR,filename,postfix);              
	sprintf(file_csv,"%s%s.%s.csv",OUTDIR,filename,postfix);   
	file_txt=malloc(sizeof(char)* (strlen(param.LAND_MODEL)+6+strlen(OUTDIR)+strlen(postfix)));	
	sprintf(file_txt,"%s%s.%s.txt",OUTDIR,param.LAND_MODEL,postfix);              

	/* 
		Draw the representation
	*/

//	printf("%f %d %d %d %d %d %f %d %d %d\n<BR>",(float)threshold, opt_log, !opt_clear, opt_ref,opt_to, GenoRef, rescale_height/100.0,R,LAND_ONLY_MUT,opt_chains);
//	outputparam(param,junk);
	

	
	draw_header();      


/*	for (i=0;i< land.nlocus;i++)
		{
		if (land.alleles[i]!=2)
		   {
		   flag=1;break;
		   }
		} */
	for (i=0;i< land.ngenotypes;i++)	
		if (land.fitness[i]==DEFAULT_FITNESS)
			{incomp=1;break;}

		
 //printf("incomp=%d flag=%d<BR>",incomp,flag);			   
	if (flag==0 && incomp==0)
	{
	
thestats=outputstats(&land, param.LAND_THRESHOLD,  param.LAND_LOG , filename,file_csv);	
char boutdligne[1024];

sprintf(boutdligne,"download CSV stat <A HREF=\"%s/%s.%s.csv\"> file</A>",OUTDIRWEB,filename,postfix);
strcat	(thestats,boutdligne);
 	if (opt_masked==1)
			{
			thestats=realloc(thestats,sizeof(char)*(strlen(thestats)+strlen("<span style=\"color: #e0141c\"><BR>Stats on whole landscape</SPAN>")+1));
			
			strcat	(thestats,"<span style=\"color: #e0141c\"><BR>Stats on whole landscape</SPAN>");
			
			}

 	
 	}
	else
/*	if (flag==1)
		{
	thestats=malloc	(sizeof(char)*strlen("No statistics avalaible when more than 2 alleles\n")+1);
	strcpy(thestats,"No statistics avalaible when more than 2 alleles\n");
	}	
	else*/
	{
	thestats=malloc	(sizeof(char)*strlen("No statistics avalaible when incomplete lanscape\n")+1);
	strcpy(thestats,"No statistics avalaible when incomplete lanscape\n");
	}

	stat_formated=malloc(sizeof(char)*(strlen(thestats)+512));
	sprintf(stat_formated,"&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<BR><BR>%s<BR><BR>",thestats);

	draw_menu (land,postfix,param,stat_formated,junk); 
		
if (opt_masked==1)
	{
	int *g;
	g=malloc(sizeof(int)*land.nlocus);

	init_landscape( &newland, land.nlocus, land.alleles);
	for (i=0;i<land.ngenotypes;i++)
		{
		int2genotype(newland,i,g);

		if (compare_geno(param.LAND_MASK,g,land.nlocus)==0)
			newland.fitness[i]=DEFAULT_FITNESS;
		else
			newland.fitness[i]=land.fitness[i];
		}
	newland.minf=	land.minf;
	newland.maxf=	land.maxf;
	p_land=&newland;	
	}
	else
	p_land=&land;



	draw_FL( 1, p_land, file_out, param.LAND_THRESHOLD, param.LAND_LOG, !param.LAND_OPTCLEAN, param.LAND_DRAWFROM,param.LAND_DRAWTOEND, param.LAND_REFERENCE, param.LAND_SCALE_H/100.0,R,param.LAND_ONLY_MUT,1,param.LAND_CHAINS,param.LAND_WINDOW_WIDTHSIZE,param.LAND_WINDOW_HEIGHTSIZE,param.LAND_COMPACT,param.LAND_FLAT,""); /* draw beautifull necklaces -- colored balls */
	

//	draw_FL( 1, p_land, file_out, threshold, opt_log, !opt_clear, opt_ref,opt_to, GenoRef, rescale_height/100.0,R,LAND_ONLY_MUT,1,opt_chains,LAND_WINDOW_WIDTHSIZE,LAND_WINDOW_HEIGHTSIZE,LAND_COMPACT,LAND_FLAT); /* draw beautifull necklaces -- colored balls */
				

//	draw_end();

if (!(flag==0 && incomp==0))
	{
	*thestats='\0';
	}	


	save_SVG( p_land,file_out, &param,R,thestats);
	save_Landscape ( p_land, 	file_txt );
	//save_SVG( p_land,file_out, threshold, opt_log, !opt_clear, opt_ref,opt_to, GenoRef, rescale_height/100.0,R,opt_chains,LAND_ONLY_MUT,LAND_COMPACT,LAND_FLAT);
	draw_lastline(param);
	menage(OUTDIR);
	
       
	free(file_out);
	free_landscape(&land);
	if (opt_masked)
		free_landscape(&newland);
	fflush (stdout);
freeparam(&param);

  exit(1);		
}








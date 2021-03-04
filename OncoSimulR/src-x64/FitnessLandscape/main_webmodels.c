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
#include "web/webparam.h"
#include "random.h"
#include "models.h"
#include "summary_statistics.h"
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
char **getMultPartData(int *nb) {
   	int i,j,k;
    int taille_chps;
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
	//f=fopen("/tmp/debug","w");
	//fprintf(f,"%d\n%s\n",content_length,cgiinput);
	//fclose (f);
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

int lengthclean(char *str)
{
	
	int i,t=strlen(str),n=0;
	
	for (i=0;i<t;i++)
		if (str[i]!=10 && str[i]!=13 && str[i]!=32)
			n++;
	
	return n;
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
		return(i);
return(-1);		
}




struct landscape ReadString( char *data, int opt_zero ){


	char letter;
	int nlocus=1,
	    l,        /* locus */
	    g;         /* genotype */
	char val[320];
	struct landscape fl;
	
	int *alleles;
	
	float fitness;
	
	int line=0,i=0;
	
	
	while( (*data == 10 || *data == 13  || *data == ' ') && *data!=0) data++ ; //1st char must be a value
	if (*data==0)
	printf("a pb occured while reading old landscape"),exit(1);
	while( (letter = (char)data[i++]) != 10 && letter!=13){ 
		if( letter == ' ' )
			nlocus++;
//			printf("%d ", data[i-1]);
	}
if (nlocus<2)
printf("a pb occured while reading old landscape nb locus ot ok"),exit(1);

	/*fprintf(stderr, "from file %s, there are %d locus\n", filename, nlocus);*/

	alleles = (int *) malloc( (size_t) nlocus * sizeof( int )  );
	if(! alleles )fprintf(stderr, "ReadFile: cannot allocate alleles, bye\n"), exit(3);

	
	for( l=0; l<nlocus; l++)
		{
		sscanf(data, "%s", val );
		*(alleles+l)=atoi(val);

		data+=strlen(val);
		while (*data==' ' && *data !='\n') data++;
		
	}


	init_landscape( &fl, nlocus, alleles);
	
	for (g=0;g<fl.ngenotypes;g++)
		if( opt_zero )
			fl.fitness[g] = 0;
		else
			fl.fitness[g] = DEFAULT_FITNESS;
		
	line = 1;
	while (*data==' ' || *data==10 || *data==13) data++; //be sure 1st char is to read
	while(*data!='\0' && *data!='-'){
		
	
		for( l=0; l<nlocus; l++)
			{
			sscanf(data, "%s", val );
		
			*(alleles+l)=atoi(val);
			data+=strlen(val);
			while (*data==' ' && *data !='\n') data++;
//			printf("%d,",atoi(val));
			}
		
			sscanf(data, "%s", val );
			fitness=atof(val);
			data+=strlen(val);
//			printf(" %f",fitness);
				
		while ((*data==' '|| *data==10 || *data==13) && *data !='\0'&& *data !='-') data++;
	
		if (*data !='\0')
		line ++;
//		printf("data=%d (%c)<BR>\n ",*data,*data);
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
		}

	/* then  find where the data begins and *so funny* this differs from one navigator to another :~/ */

	if((bout=strstr(data,"octet-stream"))==NULL)
		{
		if((bout=strstr(data,"x-fluid"))==NULL)
		printf("bad data %s <BR>\n",data),exit(1);
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



char *clean (char *data)	
{
	int i=strlen(data)-1;
	while (data[i]=='-' ||data[i]==' '||data[i]==10 ||data[i]==13) data[i--]='\0';
	i=0;
	while (*data <'L') data++ ; 
	return(data);
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





struct landscape choose_model(AllParams param,struct model_opt allopt)
//struct landscape choose_model(int which_model,int nlocus,int *nalleles,int K,int random)
{
//char *models[]={"LAND_UniRandom","LAND_Multiplicative","LAND_Kauffman_Nk","LAND_RMF","LAND_Spinglass"};

struct landscape l;
int temp;

if ( param.LAND_NBLOCI<2 )
	printf("you must have at least 2 locus <BR>"),exit(1);


init_landscape( &l, param.LAND_NBLOCI, param.ALLELS);
temp=param.LAND_LOG;
param.LAND_LOG=0;	
switch (param.LAND_MODELNBR)
	{
	case 0:
		HouseOfCards( &l ,allopt.sigma_hoc);
		break;
	
	case 1:
//		printf("model MULT.....<BR>");
		Multiplicative( &l, allopt.mu_s, allopt.sigma_s,param.LAND_VALUE_MULTSAME,allopt.DimRet,1);
		break;
	 
	case 2:
		Kaufman_NK(&l,allopt.kaufman_K, allopt.kaufman_rand);
		break;
	
	
	case 4:

		Ising(&l,allopt.mu_ising, allopt.sigma_ising, allopt.circular_ising);
		break;
		
	
	case 5:
	 	EggBox(&l,allopt.mu_eggbox,allopt.sigma_eggbox);
		break;	
	
	
	case 6:
	 	Optimum(&l,allopt.mu_prod,allopt.sigma_prod,allopt.mu_optimum,allopt.sigma_optimum);
		break;	
	
	case 3:	
	case 7:
	 	MixedModel(&l, allopt,1);
		break;	
	
	}
exp_landscape( &l );	
setminMax(&l);
param.LAND_LOG=temp;
return l;
}

void getallels(char **cgivars,AllParams param)
//void getallels(char **cgivars,int nbgeno,int *nballels,int allFields,int flag)
{
int i,a;
char suffix[32];

if (param.LAND_CHECK_SAMENBR==1)
{
a=get_value("LAND_ALL_0",cgivars,param.nbFields);

for (i=0;i<param.LAND_NBLOCI;i++)
param.ALLELS[i]=a;
}
else
for (i=0;i<param.LAND_NBLOCI;i++)
	{
	sprintf(suffix,"LAND_ALL_%d",i);
	param.ALLELS[i]=get_value(suffix,cgivars,param.nbFields);
	}
}

int getmask(char **cgivars,AllParams param)
{
int i,masked=0;
char suffix[32];

for (i=0;i<param.LAND_NBLOCI;i++)
	{
	sprintf(suffix,"LAND_MASK_%d",i);
	param.LAND_MASK[i]=get_value(suffix,cgivars,param.nbFields);
	
	if (param.LAND_MASK[i]!=-1)
		masked=1;
	}

return(masked);
}


/*void initparam(AllParams *param)
{
char *deft="0,1,1,-1,-1,0,1,1,1,0";
param->LAND_MODEL=NULL; //model name
param->LAND_NBLOCI=-1; //nbr of gentotype
param->LAND_MODELNBR=-1;
param->ALLELS=NULL; //how many allels in each genotype
param->LAND_CHECK_SAMENBR=1;	// 0 if differents 1 if same
param->LAND_LOG=0;
param->LAND_REFERENCE=-1; //nbr of left ref 
param->LAND_THRESHOLD=1; //ratio fitness
param->LAND_OPTCLEAN=0;// dsiplay neutrals
param->LAND_DRAWFROM=-1; //display lnks from
param->LAND_DRAWTOEND=-1; //display lnks from
param->LAND_MISSING=1; //1 if 0 
param->LAND_SCALE_W=1; //width scale
param->LAND_SCALE_W=1; //heigth scale
param->nbFields=0; // 
param->LAND_DEFAULT=deft;
param->LAND_RMFWEIGHT=-1;
param->LAND_PREVIOUS_FILE =NULL;
param->LAND_MULTSAME=0;
param->LAND_CHAINS=0;
param->LAND_VALUE_MULTSAME=0;
//param->LAND_SPININCOMPATIBILY=-1;
//param->LAND_RAND_SPIN=0;






param->LAND_ONLY_MUT=-1;

param->LAND_COMPACT=0;
param->LAND_COMPACT=0;
param->LAND_FLAT	=0;

}
void freeparam(AllParams *param)
{
free(param->ALLELS);
free(param->LAND_MASK);
free(param);
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

int main (int argc, char **argv)
{
	char **cgivars;
	char *thestats;

	int nbFields,masked=0;
	int ii;
	struct landscape land;
	struct landscape newland;
	struct landscape *p_land;
	int R;
	char *stat_formated;
	char *file_out=NULL;
	char *file_csv=NULL;
	char *file_txt=NULL;
	
	char postfix[1024];
	struct model_opt allopt;
	int i;//, nbmodels=5;
	
	seed_ran1( (long) time(NULL) );	
	AllParams param;
	int NBMODELS=8;
	char *models[]={"LAND_UniRandom","LAND_Multiplicative","LAND_Kauffman_Nk","LAND_RMF","LAND_Ising","LAND_EggBox","LAND_Optimum","LAND_Full_Models"};

	printf("Content-type: text/HTML\n\n");

	menage(OUTDIR);
	initparam(&param);

	cgivars=getMultPartData(&nbFields);
	param.nbFields=nbFields;
	/*leave this commented to check easily the parsing of web page*/
/*	printf("---> cgi %d\n<BR>",nbFields);
	for (ii=0;ii<nbFields*2;ii++){
	printf("%s \n",cgivars[ii]);
	if ((ii+1) %2==0)printf(" <BR>\n");
	}
	printf("end cgi\n<BR>");*/
	allopt= init_model();
	for (ii=0;ii<param.nbFields*2;ii+=2)
			{
				if 	(strcmp(cgivars[ii],"LAND_THRESHOLD")==0) 
					param.LAND_THRESHOLD=atof(cgivars[ii+1]);
											
				else	
				if 	(strcmp(cgivars[ii],"LAND_REFERENCE")==0)
					param.LAND_REFERENCE=atoi(cgivars[ii+1]);
					
				else	
				if 	(strcmp(cgivars[ii],"LAND_DRAWFROM")==0)
					param.LAND_DRAWFROM=atoi(cgivars[ii+1]);
					
				else	
				if 	(strcmp(cgivars[ii],"LAND_OPTCLEAN")==0)
					param.LAND_OPTCLEAN=1;

				else
				if 	(strcmp(cgivars[ii],"LAND_LOG")==0)
					 param.LAND_LOG=1;
					 
				else	
				if 	(strcmp(cgivars[ii],"LAND_MISSING")==0)
					 param.LAND_MISSING=1;

				else			
				if (strcmp(cgivars[ii],"LAND_SCALE_W")==0)
					param.LAND_SCALE_W=atof(cgivars[ii+1]);
											
				else	
				if (strcmp(cgivars[ii],"LAND_SCALE_H")==0)
					param.LAND_SCALE_H=atof(cgivars[ii+1]);

				else
				if (strcmp(cgivars[ii],"LAND_MODEL")==0)
					param.LAND_MODEL=clean(cgivars[ii+1]);
					
				else	
				if (strcmp(cgivars[ii],"LAND_NBLOCI" )==0)
					param.LAND_NBLOCI=atoi(cgivars[ii+1]);

				else	
				if (strcmp(cgivars[ii],"LAND_NK_PARAM")==0)
					allopt.kaufman_K=atoi(cgivars[ii+1]);

				else	
				if (strcmp(cgivars[ii],"LAND_RAND_KAUF")==0)
					allopt.kaufman_rand=1;
													
				else	
				if (strcmp(cgivars[ii],"LAND_CHECK_SAMENBR"	)==0) // si ca nest aps ds la chaine darg alors alleles differents
					param.LAND_CHECK_SAMENBR=1;
											
				else	
				if (strcmp(cgivars[ii],"LAND_PREVIOUS_FILE"	)==0) //old landscape
					{param.LAND_PREVIOUS_FILE =cgivars[ii+1];}
												
				else
				if 	(strcmp(cgivars[ii],"LAND_RMFWEIGHT"	)==0)
					param.LAND_RMFWEIGHT=atof(cgivars[ii+1]);
				else
				if 	(strcmp(cgivars[ii],"LAND_MULTSAME")==0)
					param.LAND_MULTSAME	=1;
				else
				if 	(strcmp(cgivars[ii],"LAND_DRAWTOEND")==0)
					param.LAND_DRAWTOEND=atoi(cgivars[ii+1]);	
				else
				if (strcmp(cgivars[ii],"LAND_CHAINS")==0)
					param.LAND_CHAINS=1;	
				else		
					if (strcmp(cgivars[ii],"LAND_VALUE_MULTSAME")==0)
						param.LAND_VALUE_MULTSAME=atof(cgivars[ii+1]);
				else
				 	if (strcmp(cgivars[ii],"LAND_mu_a")==0)
					{allopt.mu_s=atof(cgivars[ii+1]);}
				else	
					if (strcmp(cgivars[ii],"LAND_mu_i")==0)
					{allopt.mu_ising=atof(cgivars[ii+1]);}
				else	
					if (strcmp(cgivars[ii],"LAND_sigma_a")==0)
					{allopt.sigma_s=atof(cgivars[ii+1]);}
				else	
					if (strcmp(cgivars[ii],"LAND_sigma_h")==0)
					{allopt.sigma_hoc=atof(cgivars[ii+1]);}
				else	
					if (strcmp(cgivars[ii],"LAND_sigma_i")==0)
					{allopt.sigma_ising=atof(cgivars[ii+1]);}
				else	
					if (strcmp(cgivars[ii],"LAND_mu_e")==0)
					{allopt.mu_eggbox=atof(cgivars[ii+1]);}
				else	
					if (strcmp(cgivars[ii],"LAND_sigma_e")==0)
					{allopt.sigma_eggbox=atof(cgivars[ii+1]);}
				else	
					if (strcmp(cgivars[ii],"LAND_mu_o")==0)
					{allopt.mu_optimum=atof(cgivars[ii+1]);}
				else	
					if (strcmp(cgivars[ii],"LAND_sigma_o")==0)
					{allopt.sigma_optimum=atof(cgivars[ii+1]);}
				else	
					if (strcmp(cgivars[ii],"LAND_mu_p")==0)
					{allopt.mu_prod=atof(cgivars[ii+1]);}
				else	
					if (strcmp(cgivars[ii],"LAND_sigma_p")==0)
					{allopt.sigma_prod=atof(cgivars[ii+1]);}
					
				else	
					if (strcmp(cgivars[ii],"LAND_CIRCULAR_ISING")==0)
					{allopt.circular_ising=atof(cgivars[ii+1]);}
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
					if (strcmp(cgivars[ii],"LAND_FIX_AMOUNT")==0)
							{allopt.fix=atof(cgivars[ii+1]);}
				else	
					if (strcmp(cgivars[ii],"LAND_DIM_RET")==0)
							{allopt.DimRet=atof(cgivars[ii+1]);}

					
														
			}		
						
			param.ALLELS=malloc(sizeof(int)*param.LAND_NBLOCI);
			param.LAND_MASK=malloc(sizeof(int)*param.LAND_NBLOCI);
			getallels(cgivars,param);
//			getallels(cgivars,param.LAND_NBLOCI,nballels,nbFields,samenbr);
				
	
		masked=getmask(cgivars,param);
		
	for (i=0;i<NBMODELS;i++)
		if (strncmp(param.LAND_MODEL,models[i],strlen(models[i]))==0)
			{
			param.LAND_MODELNBR=i;break;
			}
	
		R = 7*(param.LAND_SCALE_W/100.0);
	
  		
 // 		modelname=clean (cgivars[param.LAND_MODEL]);
  			// model name et nbr???
 //  printf("---- maske %d	<BR>", masked);		
//	outputparam(param,allopt);
//	printf("mu %f sigma %f dim %f<BR>",allopt.mu_s,allopt.sigma_s,allopt.DimRet);  
	
//  printf("---- suite	<BR>");	
	if (lengthclean( param.LAND_PREVIOUS_FILE )<2) /*1st graph with this model sometimes some garbage so not0*/
		land = choose_model(param,allopt);	
	else
		{land = ReadString(param.LAND_PREVIOUS_FILE ,param.LAND_MISSING);}
//  printf("---- suite done	<BR>");	
//print_landscape(&land);
//	exit(1);//
	sprintf(postfix,"%d.%ld",getpid(),random());

	file_out=malloc(sizeof(char)* (strlen(param.LAND_MODEL)+6+strlen(OUTDIR)+strlen(postfix))); //6= 2points  and svg will be added in filename below
	file_csv=malloc(sizeof(char)* (strlen(param.LAND_MODEL)+6+strlen(OUTDIR)+strlen(postfix))); //6= 2points  and svg will be added in filename below
	file_txt=malloc(sizeof(char)* (strlen(param.LAND_MODEL)+6+strlen(OUTDIR)+strlen(postfix)));
	if(!file_out)fprintf(stderr, "main: cannot allocate file_out, bye"), exit(3);
	if(!file_csv)fprintf(stderr, "main: cannot allocate file_csv, bye"), exit(3);
	srand(time(NULL));
	//create dir where files are going to be written
	
	sprintf(file_out,"%s%s.%s.svg",OUTDIR,param.LAND_MODEL,postfix);   
	sprintf(file_csv,"%s%s.%s.csv",OUTDIR,param.LAND_MODEL,postfix);
	
	file_txt=malloc(sizeof(char)* (strlen(param.LAND_MODEL)+6+strlen(OUTDIR)+strlen(postfix)));	
	sprintf(file_txt,"%s%s.%s.txt",OUTDIR,param.LAND_MODEL,postfix);              
	param.LAND_OUTFILE=malloc (sizeof(char )*(strlen(file_out)  +1));
	strcpy(param.LAND_OUTFILE,file_out);
	

	/* 
		Draw the representation
	*/
	draw_header(1000,800); 
	

	//thestats= outputstats(&land,1);
	/*if (param.LAND_CHECK_SAMENBR==1 && param.ALLELS[0]==2 && (land.minf!=land.maxf))
		{
		thestats=outputstats(&land, param.LAND_THRESHOLD, param.LAND_LOG , models[param.LAND_MODELNBR]);
		if (masked==1)
			{
			thestats=realloc(thestats,sizeof(char)*(strlen(thestats)+strlen("<span style=\"color: #e0141c\"><BR>Stats on whole landscape</SPAN>")+1));
			strcat	(thestats,"<span style=\"color: #e0141c\"><BR>Stats on whole landscape</SPAN>");
			}
		}
	else
		thestats=malloc	(sizeof(char)*strlen("For now, no statistics avalaible for allels >2 or flat landscapes (same fitness everywhere)\n")+1);
		strcpy(thestats,"For now, no statistics avalaible <BR>for allels >2 or flat landscapes\n");
*/
if (land.minf!=land.maxf)	
	{
		thestats=outputstats(&land, param.LAND_THRESHOLD, param.LAND_LOG , models[param.LAND_MODELNBR],file_csv);
//		strcat	(thestats,"download CSV stat <A HREF=\"%s/%s.csv\"> file</A><BR>",OUTDIRWEB,param.LAND_MODEL,postfix);

		char boutdligne[1024];

		sprintf(boutdligne,"download CSV stat <A HREF=\"%s/%s.%s.csv\"> file</A>",OUTDIRWEB,param.LAND_MODEL,postfix);
		strcat	(thestats,boutdligne);



		if (masked==1)
			{
			thestats=realloc(thestats,sizeof(char)*(strlen(thestats)+strlen("<span style=\"color: #e0141c\"><BR>Stats on whole landscape</SPAN>")+1));
			strcat	(thestats,"<span style=\"color: #e0141c\"><BR>Stats on whole landscape</SPAN>");
			}
		}
else
{
	thestats=malloc	(sizeof(char)*strlen("No statistics avalaible for  flat landscapes (same fitness everywhere)\n")+1);
		strcpy(thestats,"No statistics avalaible for flat landscapes (same fitness everywhere)\n");
}	
			
//	draw_menu (land,cgivars,nbFields,ifile+1,postfix,param,thestats,allopt); 
	stat_formated=malloc(sizeof(char)*(strlen(thestats)+512));
	sprintf(stat_formated,"<BR>%s<BR><BR>",thestats);
	draw_menu (land,postfix,param,stat_formated,allopt); 
	
	//if mask is on then a newlandscape has to born
if (masked==1)
	{
	int *g;
	g=malloc(sizeof(int)*param.LAND_NBLOCI);

	init_landscape( &newland, param.LAND_NBLOCI, param.ALLELS);
	for (i=0;i<land.ngenotypes;i++)
		{
		int2genotype(newland,i,g);

		if (compare_geno( param.LAND_MASK,g,param.LAND_NBLOCI)==0)
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
		
//	draw_end();

	draw_FL( 1, p_land, file_out, param.LAND_THRESHOLD, param.LAND_LOG, !param.LAND_OPTCLEAN, param.LAND_DRAWFROM,param.LAND_DRAWTOEND, param.LAND_REFERENCE, param.LAND_SCALE_H/100.0,(int)R,param.LAND_ONLY_MUT,1,param.LAND_CHAINS,param.LAND_WINDOW_WIDTHSIZE,param.LAND_WINDOW_HEIGHTSIZE,param.LAND_COMPACT,param.LAND_FLAT,""); /* draw beautifull necklaces -- colored balls */
	if (land.minf==land.maxf)
		*thestats='\0';
//	save_SVG( p_land,file_out, param.LAND_THRESHOLD, param.LAND_LOG, param.LAND_OPTCLEAN, param.LAND_DRAWFROM,param.LAND_DRAWTOEND, param.LAND_REFERENCE, param.LAND_SCALE_H/100.0,R,param.LAND_CHAINS,param.LAND_ONLY_MUT,param.LAND_COMPACT,param.LAND_FLAT);
		save_SVG( p_land,file_out, &param,R,thestats);
	save_Landscape ( p_land, 	file_txt );
//outputparam(param,allopt);
	draw_lastline(param);
	
	
	
	
       
	free(file_out);

//print_landscape(p_land);
freeparam(&param);
free_landscape(&land);
if (masked)
free_landscape(&newland);


  exit(1);		
}








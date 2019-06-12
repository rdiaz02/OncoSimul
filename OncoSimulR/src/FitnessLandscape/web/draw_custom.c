#include <stdio.h>
#include <stdlib.h>
#include  <string.h>
#include <unistd.h>
#include <float.h>

#include "landscape.h"
#include "genotypes.h"
#include "drawings.h"
#include "random.h"
#include "models.h"
#include "summary_statistics.h"
#include "web/webparam.h"
#include "web/web_common.h"
#include "web/server.h"
//on my mac




/*****************************/
int main(int argc,char **argv)
{
	FILE *f;
	char *data;
	int nlog,i,incomp=0,flag=0;;
	char filename[1032],*thestats,*p;
	struct landscape land;
	char where[132],*file_out,*file_csv;
	
	AllParams params;
	struct model_opt junk;
	int R=5;	
		
	printf("%s%c%c\n","Content-Type:text/html;charset=iso-8859-1",13,10);
	data = getenv("QUERY_STRING");
	if(data == NULL)
  		{printf("<P>Error! Error in passing data from form to script %s",data);exit(1);}
	else 
		if (strchr(data,'&')==NULL)
  			{printf("<P>Error! Invalid data %s ",data);exit(1);}
  	if ((strstr(data,"LAND_FILE")==NULL) ||	(strstr(data,"LAND_LOG")==NULL	))
   		{printf("<P>Error! Error in passing data from form to script %s",data);exit(1);}
	
	
	initparam(&params);
  	p=strstr(data,"LAND_LOG=");
  	if (p==NULL)  	
  	printf("<P>Error! Error in passing form (log value not found)"),exit(1);
  	p+=strlen("LAND_LOG=");

  	
  	nlog=*p -'0';
  	sprintf(filename,"%s",LAND_FILE_PLACE);
  	p=strstr(data,"LAND_FILE=");
  	 	if (p==NULL)  	
  	printf("<P>Error! Error in passing form (filename value not found)"),exit(1);
  	p+=strlen("LAND_FILE=");
  	
  	params.LAND_MODEL=malloc(sizeof(char)*(strlen(p)+1));    //copy file name instaed of model name    
	strcpy(params.LAND_MODEL,p);

  	strcat(filename,p);

  	sprintf(where,"%d.%ld",getpid(),random());	

 
  	
  	f=fopen(filename,"r");//test the file
  	if (f==NULL)
  	printf("File not found . Please mail the author(s)s name(s) to sophieb(at)abi.snv.jussieu.fr <BR>"),exit(0);
  	fclose(f);
 	land = ReadFile(filename,0);
 /*	for (i=0;i<land.nlocus-1;i++)
 		{
 		if (land.alleles[i]!=land.alleles[i+1])
 			{flag=1;break;}
 		}*/
 	setminMax(&land);
 	
	params.LAND_NBLOCI=land.nlocus;
	
	params.ALLELS=malloc(sizeof(int)* land.nlocus);
	
	for (i=0;i< land.nlocus;i++)	
		params.ALLELS[i]=land.alleles[i];
		
	params.LAND_MASK=malloc(sizeof(int)*land.nlocus);
	for (i=0;i< land.nlocus;i++)	
		params.LAND_MASK[i]=-1;	
		
 	for (i=0;i< land.ngenotypes;i++)	
			if (land.fitness[i]==DEFAULT_FITNESS)
				{incomp=1;break;}

 	sprintf(where,"%d.%ld",getpid(),random());

	file_out=malloc(sizeof(char)* (strlen(params.LAND_MODEL)+16+strlen(OUTDIR)+strlen(where))); //6= 2points  and svg will be added in filename below
	file_csv=malloc(sizeof(char)* (strlen(params.LAND_MODEL)+16+strlen(OUTDIR)+strlen(where))); //6= 2points  and csv will be added in filename below
	
	
	if(!file_out)fprintf(stderr, "main: cannot allocate file_out, bye"), exit(3);
	sprintf(file_out,"%s%s.%s.svg",OUTDIR,params.LAND_MODEL,where);  
	sprintf(file_csv,"%s%s.%s.csv",OUTDIR,params.LAND_MODEL,where);
	
	draw_header();
	if (incomp==1)
		{
		thestats=malloc	(sizeof(char)*strlen("No statistics avalaible when incomplete lanscape\n")+1);
		strcpy(thestats,"No statistics avalaible when incomplete lanscape\n");
		}
	else
	/*	if (flag==1)
			{
			thestats=malloc	(sizeof(char)*strlen("No statistics avalaible when more than 2 alleles\n")+1);
			strcpy(thestats,"No statistics avalaible when more than 2 alleles\n");
			}
		else*/
			{
			thestats=NULL;
			thestats=outputstats(&land, 1, nlog , filename,file_csv);	
			
			}
		
	
 	
 
	 params.LAND_THRESHOLD=1;
	 params.LAND_LOG=nlog;
	 params.LAND_OPTCLEAN=1;
	 params.LAND_DRAWFROM=-1;
	 params.LAND_DRAWTOEND=-1;
	 params.LAND_REFERENCE=0;
	params.LAND_ONLY_MUT=-1;
	params.LAND_CHAINS=0;
	params.LAND_COMPACT=0;
	params.LAND_FLAT=0;
	
	draw_menu (land,where,params,thestats,junk); 	
	if ((incomp==1) || (flag==1))
		*thestats='\0';
	draw_FL( 1, &land, file_out, 1,nlog, 1, -1,-1, 0, 1.0 ,R,  -1,  1,0,  -1, -1,0,0,NULL); /* draw beautifull necklaces -- colored balls */

//	draw_FL( 0, fl, file_out, param->LAND_THRESHOLD, param->LAND_LOG, !param->LAND_OPTCLEAN, param->LAND_DRAWFROM,param->LAND_DRAWTOEND, param->LAND_REFERENCE, param->LAND_SCALE_H/100.0,R, param->LAND_ONLY_MUT,1,param->LAND_CHAINS,1000,1000,param->LAND_COMPACT,param->LAND_FLAT);         


	save_SVG(&land,file_out, &params,R,thestats);
		
	draw_lastline(params);
	
	
 	free_landscape(&land);
 	freeparam(&params);
 	return(0);
}

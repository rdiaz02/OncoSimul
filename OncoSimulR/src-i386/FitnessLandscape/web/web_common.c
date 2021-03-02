#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <float.h>



#include "landscape.h"
#include "models.h"
#include "webparam.h"
#include "web_common.h"
#include "drawings.h"
#include "server.h"

//on my mac




void setminMax(struct landscape *l)
{
int g;
float min=FLT_MAX;
float max=FLT_MIN;
for(g=0; g<l->ngenotypes; g++ )
	{
	if (l->fitness[g]!=DEFAULT_FITNESS)
		{
		if (l->fitness[g]>max)
			max=l->fitness[g];
		if (l->fitness[g]<min)
			min=l->fitness[g];
		
		}
	}	
l->minf=min;
l->maxf=max;	
}


 void draw_header()
{
	printf("<!DOCTYPE html>\n<html>\n");
	printf("<head>\n");
//	printf("<meta http-equiv=\"content-type\" content=\"application/xhtml+xml; charset=utf-8\" />\n");
	printf("<title>MAGELLAN</title>\n");
	printf("<meta content=\"landscape manipulation\" name=\"keywords\" />");
    printf("	<style type=\"text/css\" media=\"all\">\n");
    printf("          @import \"%smagellan.css\";\n",MAINDIR);
 	printf("	</style>\n");
 	printf("<script type=\"text/javascript\" src=\"%slandmodel_js.js\">",JSDIR );
	
	printf("</script>");

	printf(  "</head><body>\n");
/*
	printf(" <div id=\"haut2\">\n");
	printf("	<table>\n");
	printf("	<tr>\n");
	printf("	<td>\n");
	printf("	<A href=\"%s%s\"><H5>MAGELLAN main </H5></A>\n",INITDIRWEB,HTML_MAIN);

	printf("	</td><td align=\"right\"></td></tr>\n");
	printf("	</table>\n");
	printf("</div >\n");*/
}

void draw_end()
{
printf(  "</html>\n");
printf(  "<!-- end -->\n");
}

/*
IMPORTRANT check if the number of the model is always uptodate
LAND_UniRandom=0,LAND_Multiplicative=1,LAND_Kauffman_Nk=2,LAND_RMF=3,LAND_Ising=4,LAND_EggBox=5,LAND_Optimum=6,LAND_FULL=7;
*/
void write_html_model( AllParams  param,struct model_opt allopt)
{
	
	switch (param.LAND_MODELNBR)
	{
	case 0: //(m=="LAND_UniRandom")
			printf(	"<table><TR><th>Name</TH><th>mean</TH><th>stdev</TH></TR>\n");
			printf(	"<TR><td>HOC</Td><td>0</Td><td><input type=\"text\" value=\"%f\" size=3 id=\"LAND_sigma_h\" name=\"LAND_sigma_h\" readonly></Td>\n",allopt.sigma_hoc);
			
			printf(	"</TR>\n");
			printf(	"</Table>\n");
			
	break;
	
	case 1: 	//LAND_Multiplicative
		
			printf("<table><TR><th>Name</TH><th>mean</TH><th>stdev</TH><TH></TH><TH></TH></TR><TR>\n");
			printf(	"<TR><td>s:</Td><td><input type=\"text\" value=\"%f\" size=3 id=\"LAND_mu_a\"  name=\"LAND_mu_a\" readonly ></Td><td><input type=\"text\" value=\"%f\" size=3 id=\"LAND_sigma_a\" name=\"LAND_sigma_a\" readonly  ></Td>\n",allopt.mu_s,allopt.sigma_s);
			printf(	"<TD>DimRet:</TD><TD><Input type=\"text\" value=\"%f\" size=3 id=\"LAND_DIM_RET\" name=\"LAND_DIM_RET\" readonly  STYLE=\"color: lightgrey;\" ></TD>\n",allopt.DimRet);
			printf(	"</TR>\n");
			printf(	"</Table>\n");
		
	break;
	
	case 2://(m=="LAND_Kauffman_Nk")
		
			printf(	"<table><TR><th></TH><th>loci</TH><th>random</TH></TR>\n");
			printf(	"<TR><td>Kauf</Td><td><input type=\"text\" value=\"%d\" size=3 id=\"LAND_NK_PARAM\" name=\"LAND_NK_PARAM\" readonly >",allopt.kaufman_K);
			printf(	"</Td><td><input type= \"checkbox\" id=\"LAND_RAND_KAUF\" name=\"LAND_RAND_KAUF\" readonly ");
			if (allopt.kaufman_rand==1) printf("checked");
			printf(	"	></Td>\n");
			printf(	"</TR>");
			printf(	"</table>\n");
	break;

	case 3: //(m=="LAND_RMF")
		
			printf(	"<table><TR><th>Name</TH><th>mean</TH><th>stdev</TH><th>Dim Ret</TH></TR>\n");
			printf(	"<TR><td>HOC</Td><TD></Td><td><input type=\"text\" value=\"%f\" size=3 id=\"LAND_sigma_h\" name=\"LAND_sigma_h\" ></Td><TD></TD>\n",allopt.sigma_hoc);
			printf(	"</TR>\n");
			printf(	"<TR><td>Mult</Td><td><input type=\"text\" value=\"%f\" size=3 id=\"LAND_mu_a\" name=\"LAND_mu_a\" ></Td><td><input type=\"text\" value=\"%f\" size=3 id=\"LAND_sigma_a\" name=\"LAND_sigma_a\"></Td>\n",allopt.mu_s,allopt.sigma_s);
			printf(	"<TD><Input type=\"text\" value=\"%f\" size=3 id=\"LAND_DIM_RET\" name=\"LAND_DIM_RET\"></TD>\n",allopt.DimRet);
			printf(	"</TR>\n");
			printf(	"</Table>\n");
	break;

	case 4: //(m=="LAND_Ising")
			
			printf(	"<table><TR><th>Name</TH><th>mean</TH><th>stdev</TH><th>circular<TH></TR>\n");
			printf(	"<TR><td>Ising</Td><td><input type=\"text\" value=\"%f\" size=3 id=\"LAND_mu_i\" name=\"LAND_mu_i\" readonly></Td><td><input type=\"text\" value=\"%f\" size=3 id=\"LAND_sigma_i\" name=\"LAND_sigma_i\" readonly></Td>\n",allopt.mu_ising,allopt.sigma_ising);
			
			printf(	"<TD><input type= \"checkbox\" id=\"LAND_CIRCULAR_ISING\" name=\"LAND_CIRCULAR_ISING\"  readonly");
			if (allopt.circular_ising)
		 		printf("checked");
			printf(	"></TD>	\n");
			
			
			printf(	"</TR>\n");
		   printf(	"</Table>\n");

	break;
		
	case 5: //(m=="LAND_EggBox")
			
			printf(	"<table><TR><th>Name</TH><th>mean</TH><th>stdev</TH><th></TH></TR>\n");
			printf(	"<TR><td>EggBox</Td><td><input type=\"text\" value=\"%f\" size=3 id=\"LAND_mu_e\"  name=\"LAND_mu_e\"  readonly></Td><td><input type=\"text\" value=\"%f\" size=3 id=\"LAND_sigma_e\"  name=\"LAND_sigma_e\"  readonly></Td>\n",allopt.mu_eggbox,allopt.sigma_eggbox);
			printf(	"<TD></TD>\n");	

			printf(	"</TR>");
			printf(	"</table>\n");
			
			break;
	case 6: //(m=="LAND_Optimum")
		
		printf("<TABLE><TR><td></Td><TD colspan=2><b>Production</b></TD><TD colspan=2><b>Fitness Funct.</b></TD></TR>");
			printf("<TR><td></Td><td>mean</Td><td>stdev</Td><td>mean</Td><td>stdev </Td></TR>");
					printf(	"<TR><td></Td><td><input type=\"text\" value=\"%f\" size=3 id=\"LAND_mu_o\" name=\"LAND_mu_o\" readonly></TD><TD><input type=\"text\" value=\"%f\" size=3 id=\"LAND_mu_p\"  name=\"LAND_mu_p\"readonly></Td>\n",allopt.mu_optimum,allopt.mu_prod);
			printf(	"<td><input type=\"text\" value=\"%f\" size=3 id=\"LAND_sigma_o\"  name=\"LAND_sigma_o\" readonly></Td></Td><td><input type=\"text\" value=\"%f\" size=3 id=\"LAND_sigma_p\" name=\"LAND_sigma_p\" readonly></Td>\n",allopt.sigma_optimum,allopt.sigma_prod);
			printf(	"</table>\n");
	break;
	
	case 7: //(m=="LAND_Full_Modelss")
			printf("Basal Fitness: <Input type=\"text\" value=\"%f\" size=3 id=\"LAND_FIX_AMOUNT\" readonly name=\"LAND_FIX_AMOUNT\">",allopt.fix);
			printf(	"<table><TR><th>Name</TH><th>mean</TH><th>stdev</TH><th>various</TH></TR>\n");
			printf(	"<TR><td>HOC</Td><td>0</Td><td><input type=\"text\" value=\"%f\" size=3 id=\"LAND_sigma_h\"  name=\"LAND_sigma_h\" readonly></Td><td></td>\n",allopt.sigma_hoc);
			printf(	"</TR>\n");
			printf("<tr bgcolor=\"#DFDFDF\"><td></td> <td></td> <td></td> <td></td><td></td></tr>\n");
			printf(	"<tr><td></td><td></td><td></td><td></td><td></td><td></td></tr>\n");
			printf(	"<TR><td>s:</Td><td><input type=\"text\" value=\"%f\" size=3 id=\"LAND_mu_a\" name=\"LAND_mu_a\" readonly></Td><td><input type=\"text\" value=\"%f\" size=3 id=\"LAND_sigma_a\" name=\"LAND_sigma_a\" readonly></Td>\n",allopt.mu_s,allopt.sigma_s);
			printf(	"<TD>DimRet:</TD><TD><Input type=\"text\" value=\"%f\" size=3 id=\"LAND_DIM_RET\"  name=\"LAND_DIM_RET\" readonly STYLE=\"color: lightgrey;\"></TD>\n",allopt.DimRet);			
			printf(	"</TR>\n");
			printf("<tr bgcolor=\"#DFDFDF\"><td></td> <td></td> <td></td> <td></td><td></td></tr>\n");

			printf(	"<tr><td></td><td></td><td></td><td>circular</td></tr><td></td>\n");
			printf(	"<TR><td>Ising</Td><td><input type=\"text\" value=\"%f\" size=3 id=\"LAND_mu_i\" name=\"LAND_mu_i\"></Td><td><input type=\"text\" value=\"%f\" size=3 id=\"LAND_sigma_i\" name=\"LAND_sigma_i\"></Td>\n",allopt.mu_ising,allopt.sigma_ising);
			printf(	"<TD><input type= \"checkbox\" id=\"LAND_CIRCULAR_ISING\" name=\"LAND_CIRCULAR_ISING\"  readonly");
			if (allopt.circular_ising)
		 		printf("checked");
			printf(	"></TD>	\n");
			printf("<tr bgcolor=\"#DFDFDF\"><td></td> <td></td> <td></td> <td></td><td></td></tr>\n");

			printf(	"<tr><td></td><td></td><td></td><td>min</td><td></td></tr>\n");
			printf(	"<TR><td>EggBox</Td><td><input type=\"text\" value=\"%f\" size=3 id=\"LAND_mu_e\"  name=\"LAND_mu_e\" readonly></Td><td><input type=\"text\" value=\"%f\" size=3 id=\"LAND_sigma_e\" name=\"LAND_sigma_e\" readonly></Td>\n",allopt.mu_eggbox,allopt.sigma_eggbox);
			printf(	"<TD></TD>\n");	
			printf(	"</TR>");
			printf("<tr bgcolor=\"#DFDFDF\"><td></td> <td></td> <td></td> <td></td><td></td></tr>\n");

		printf("<TR><td></Td><TD colspan=2><b>Production</b></TD><TD colspan=2><b>Fitness Funct.</b></TD></TR>");
			printf("<TR><td></Td><td>mean</Td><td>stdev</Td><td>mean</Td><td>stdev </Td></TR>");
					printf(	"<TR><td></Td><td><input type=\"text\" value=\"%f\" size=3 id=\"LAND_mu_o\" name=\"LAND_mu_o\" readonly></TD><TD><input type=\"text\" value=\"%f\" size=3 id=\"LAND_mu_p\"  name=\"LAND_mu_p\"readonly></Td>\n",allopt.mu_optimum,allopt.mu_prod);
			printf(	"<td><input type=\"text\" value=\"%f\" size=3 id=\"LAND_sigma_o\"  name=\"LAND_sigma_o\" readonly></Td></Td><td><input type=\"text\" value=\"%f\" size=3 id=\"LAND_sigma_p\" name=\"LAND_sigma_p\" readonly></Td>\n",allopt.sigma_optimum,allopt.sigma_prod);
	


					printf("<tr bgcolor=\"#DFDFDF\"><td></td> <td></td> <td></td> <td></td><td></td></tr>\n");
			
			
			printf("<tr bgcolor=\"#777777\"><td></td> <td></td> <td></td> <td></td><td></td></tr>");	
			printf("<tr><TD></TD><td><b>fitness</b></Td><td></Td><td></Td><td></Td></TR>");			
			printf("<tr><TD>Minimum</TD><td> <Input type=\"text\" value=\"0\" size=3 id=\"LAND_FIX_AMOUNT\" name=\"LAND_FIX_AMOUNT\"></td><td></Td></TR>");
			
			
		printf("<tr bgcolor=\"#777777\"><td></td> <td></td> <td></td> <td></td><td></td></tr>\n");

			printf(	"<tr><td></td><td><b>loci</b></td><td><b>random</b></td><td></td><td></td></tr>\n");

				printf(	"<TR><td>Kauf</Td><td><input type=\"text\" value=\"%d\" size=3 id=\"LAND_NK_PARAM\" name=\"LAND_NK_PARAM\" readonly >",allopt.kaufman_K);
				
			printf(	"</Td><td><input type= \"checkbox\" id=\"LAND_RAND_KAUF\"  name=\"LAND_RAND_KAUF\"readonly ");
			if (allopt.kaufman_rand==1) printf("checked");
			printf(	"	></Td>\n");
		printf("<tr bgcolor=\"#777777\"><td></td> <td></td> <td></td> <td></td><td></td></tr>\n");
		
			printf(	"<td></td></TR>");
			
			printf(	"</table>\n");

			break;
		default:	
			printf(	"**PB**");
		}	


	
}

void draw_lastline(AllParams param)
{
	if (param.LAND_MODELNBR==-1)
		printf("<form name=\"land_form2\" method=\"post\" action=\"%s\" enctype=\"multipart/form-data\">\n",PROGFILE);
	else
		{
		printf("<form name=\"land_form2\" method=\"post\" action=\"%s\" enctype=\"multipart/form-data\">\n",PROGMODELE);
		printf(	"<input type=\"Button\" value= \"Generate again\" title=\"New realisation\" onClick=\"reset_previous();land_form.submit();\" />	");
		}
	printf(	"	<input type=\"Button\" value= \"Draw\" onClick=land_form.submit(); />");		
	printf("</form>\n");
	printf(  "</html>\n");
printf(  "<!-- end -->\n");

}

void draw_menu (struct landscape land,  char *where, AllParams param, char *st,struct model_opt allopt)
{	

	char *models[]={"LAND_UniRandom","LAND_Multiplicative","LAND_Kauffman_Nk","LAND_RMF","LAND_Ising","LAND_EggBox","LAND_Optimum","LAND_Full_Models"};
	int i;
		
	printf("<UL id=\"menu\">\n");
 	printf("<LI> <h3><A href=\"#\"  onclick =\"cachedecache('menu_land');myhide('menu_save');myhide('view_land');\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Landscape&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</a></h3>\n");

	printf("<div id=\"menu_land\" class=\"cache\">\n");
//	printf(	"model no%d %s<br>",param.LAND_MODELNBR,param.LAND_MODEL);exit(1);
	if (param.LAND_MODELNBR==-1)//file
	{
		printf("<form name=\"land_form\" method=\"post\" action=\"%s\" enctype=\"multipart/form-data\">\n",PROGFILE);

	printf(	"<P style=\"text-align:center;\"><input type=\"text\" name=\"LAND_FILE\" value=\"%s\" style=\"background-color:#EEEEEE;text-align:center;\" readonly size=%d></p><p>\n",
	param.LAND_MODEL,(int)strlen(param.LAND_MODEL));
	printf(	"<input type=\"Button\" value= \"New\" onClick=\"location.href='%s%s';\">\n",INITDIRWEB,HTML_MAIN);
	

	}
	else
	{
		printf("<form name=\"land_form\" method=\"post\" action=\"%s\" enctype=\"multipart/form-data\">\n",PROGMODELE);

		printf(	"<select id =\"LAND_MODEL\" name=\"LAND_MODEL\"  >\n");
		if (param.LAND_MODELNBR==0)
		printf(	"<option value=\"%s\" selected> House of cards</option>",models[param.LAND_MODELNBR]);
		else	
		printf(	"<option value=\"%s\" selected> %s</option>",models[param.LAND_MODELNBR],strchr(models[param.LAND_MODELNBR],'_')+1);
		printf(	"</select>\n");
		write_html_model(param,allopt);	
	

		printf(	"<p id=\"nbgeno\">\n");
		printf("Number of loci:");
//		n= get_value("LAND_NBLOCI",params,nb_couples);
		printf("<input type=\"text\" id=\"LAND_NBLOCI\" name=\"LAND_NBLOCI\" size=\"1\" value=\"%d\"  readonly >",land.nlocus);
		printf(	"</p>");
		
		printf(	"<p id=\"tablespace\">\n");
		printf("Enter your number of alleles <table><tbody><tr>");
		for (i=0;i<land.nlocus;i++)
			{
			char name[32];
			sprintf(name,"LAND_ALL_%d",i);
			printf("<td><input type=\"text\" id=\"%s\" name=\"%s\" size=\"1\"  value=\"%d\"  readonly></td>\n",
					name,name,land.alleles[i]);
			}
		printf("</table>");
//do not write this if same nbr
		if (param.LAND_CHECK_SAMENBR==1)
		printf(" <input type= \"hidden\" id=\"LAND_CHECK_SAMENBR\" name=\"LAND_CHECK_SAMENBR\"   value=\"0\" >");
			
		printf(	"	<input type=\"Button\" value= \"Generate again\" title=\"New realisation\" onClick=\"reset_previous();land_form.submit();\" />	");
		printf(	"<input type=\"Button\" value= \"New model\" title=\"New model\" onClick=\"location.href='%s%s';\">\n",INITDIRWEB,HTML_MODEL);
	}
	printf(" </div>\n");
	printf(" </li>\n");
	
 	printf("<li><h3><A href=\"#\"  onclick =\"myhide('menu_land');myhide('menu_save');cachedecache('view_land');\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Visual&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</a></h3>\n");

	printf("<div id=\"view_land\" class=\"cache\">\n");
	printf(	"	<br><p>\n");
	
	printf(	"<p>Reference on left: <input type=\"text\" name=\"LAND_REFERENCE\" size=\"2\" value=\"%d\"> (click on graph)</p><p>\n",param.LAND_REFERENCE);
	if (param.LAND_COMPACT==1)
		printf("<p>Compact view:<input type=\"checkbox\" name=\"LAND_COMPACT\" id=\"LAND_COMPACT\" checked><p>\n");
	else
		printf("<p>Compact view:<input type=\"checkbox\" name=\"LAND_COMPACT\" id=\"LAND_COMPACT\"><p>\n");
	if (param.LAND_FLAT==1)
		printf("<p>Flat view:<input type=\"checkbox\" name=\"LAND_FLAT\" id=\"LAND_FLAT\" checked><p>\n");
	else
		printf("<p>Flat view:<input type=\"checkbox\" name=\"LAND_FLAT\" id=\"LAND_FLAT\"><p>\n");
	
	printf(	"<HR>");
	
	printf(	"Scale:<input type=\"text\" name=\"LAND_SCALE_W\" size=\"6\" value=\"%.1f\">%%<p>\n",param.LAND_SCALE_W);
	printf(	"Heigth Scale:<input type=\"text\" name=\"LAND_SCALE_H\" size=\"6\" value=\"%.1f\">%%<p>\n",param.LAND_SCALE_H);
	if (param.LAND_LOG==1)
			printf(	"Log Fitness:<input type=\"checkbox\" name=\"LAND_LOG\" checked ><p>\n");
	else
			printf(	"Log Fitness:<input type=\"checkbox\" name=\"LAND_LOG\"  onclick =\"checklog(%f);\"><p>\n",land.minf);
	printf(	"<HR>");
	printf(	"<p><center><b>Neutral mutations </center></b>");

	printf(	"<p>Ratio:<input type=\"text\" name=\"LAND_THRESHOLD\" size=\"3\" value=\"%.1f\"  />",param.LAND_THRESHOLD);
	if (param.LAND_OPTCLEAN==1)
			printf(	"Display:<input type=\"checkbox\" name=\"LAND_OPTCLEAN\" checked> \n");
	else
			printf(	"Display:<input type=\"checkbox\" name=\"LAND_OPTCLEAN\"> <p>\n");
	printf(	"<HR>");
	printf(	"<p><center><b>Select path </center></b>");

	printf(	"<p>&nbsp;&nbsp;From: <input type=\"text\" title=\"default is -1\" name=\"LAND_DRAWFROM\" size=\"2\" value=\"%d\"> (alt-click on graph)</p>\n",param.LAND_DRAWFROM);	
	printf(	"<p>&nbsp;&nbsp;To:&nbsp;&nbsp; <input type=\"text\" title=\"default is -1\" name=\"LAND_DRAWTOEND\" size=\"2\" value=\"%d\"> (shift-click on graph)</p>\n",param.LAND_DRAWTOEND);	
	printf("<P>\n");
	if (param.LAND_CHAINS==1)
	printf("Draw only chains <input type=\"checkbox\" name=\"LAND_CHAINS\" id=\"LAND_CHAINS\" checked><p>\n");
	else
	printf("Draw only chains <input type=\"checkbox\" name=\"LAND_CHAINS\" id=\"LAND_CHAINS\" ><p>\n");
	printf("<p>Display only locus<input type=\"text\" name=\"LAND_ONLY_MUT\" size=\"1\" value=\"%d\" title=\"default is -1\"></p>\n",param.LAND_ONLY_MUT);
	printf(	"<HR>");
	printf(	"<p> <center><b>Sub landscape </center></b>");	
	printf("<p>Fix allele at locus<TABLE>");
	// ici ajouer un script qui verifie que si on change une case ca reste un chiffre possible (en fn  du tableau des alleles)
	for (i=0;i<land.nlocus;i++)
		printf("<th>%d</th>",i+1);
	printf("<tr>");	
	
	for (i=0;i<land.nlocus;i++)
			{
			char name[32];
			sprintf(name,"LAND_MASK_%d",i);
			
			printf("<td><input type=\"text\" id=\"%s\" name=\"%s\" size=\"1\"  value=\"%d\">  </td>\n",
					name,name,(int)param.LAND_MASK[i]);
			}
			printf("</TR></table>");
		printf(	"<HR>");
			

	printf(	"<input type=\"Button\" value= \"Clear\" onClick=\"resetValues(%d)\">\n",land.nlocus);

	printf(	"<input type =\"hidden\" name=\"LAND_PREVIOUS_FILE\" value=\"");
	output_landscape( &land);

	printf("\" >\n");
	
	printf(	"	<input type=\"Button\" value= \"Draw\" onClick=land_form.submit(); />	");


	printf(	"</form>\n");
	printf("</div>\n");
	printf("</LI>\n");

		printf(	"<LI ><h3><A Href=\"#\"   onclick =\"myhide('menu_land');myhide('view_land');cachedecache('menu_save');\">&nbsp;&nbsp;Save&nbsp;&nbsp;</a></h3>\n");
	printf(	"	\n");
	printf("		<div id=\"menu_save\" class =\"cache\">\n");
	printf(" Click selected format:<BR>");
	if(param.LAND_MODELNBR!=-1)//model)
	{

	printf("<br><A Href=\"%s%s.%s.svg\"  download=\"%s.svg\">Save SVG</a>",OUTDIRWEB,models[param.LAND_MODELNBR],where,models[param.LAND_MODELNBR]);
	printf("<br><A Href=\"%s%s.%s.pdf\"  download=\"%s.pdf\">Save PDF</a>",OUTDIRWEB,models[param.LAND_MODELNBR],where,models[param.LAND_MODELNBR]);
	printf("<br><A Href=\"%s%s.%s.txt\"  download=\"%s.txt\">Save TXT</a>",OUTDIRWEB,models[param.LAND_MODELNBR],where,models[param.LAND_MODELNBR]);
	
	}
	else
	{
//	printf("<br><A Href=\"%s%s.%s.txt\"  download=\"%s.txt\">Save TXT*</a>",OUTDIRWEB,models[param.LAND_MODELNBR],where,models[param.LAND_MODELNBR]);

	printf("<br><A Href=\"%s%s.%s.svg\"  download=\"%s.svg\">Save SVG</a>",OUTDIRWEB,param.LAND_MODEL,where,param.LAND_MODEL);
	printf("<br><A Href=\"%s%s.%s.pdf\"  download=\"%s.pdf\">Save PDF</a>",OUTDIRWEB,param.LAND_MODEL,where,param.LAND_MODEL);
	
	}
	printf("		<BR></div>\n");
	printf(	"	\n");

	printf("</LI>\n");



	
/*	printf(	"<LI ><h3><A href=\"#\" onclick =\"cachedecache(\'LAND_STAT\');myhide('menu_land');myhide('view_land');myhide('menu_save');\">&nbsp;&nbsp;&nbsp;&nbsp;Stats&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp</a></h3>\n");
	printf("<p id=\"LAND_STAT\" class =\"cache\">%s",st); 	printf(	"<br><br><br><br>\n");
	printf(	"</P></LI >");*/
	printf(	"<LI ><h3><A Href=\"%s%s\" target=\"_blank\">Help</a></h3></LI >\n",INITDIRWEB,HTML_HELP);
	printf(	"<LI ><h3><A Href=\"%s%s\" >Home</a></h3></LI >\n",INITDIRWEB,HTML_MAIN);
	
		printf(	"</UL>\n");
	printf(	"<br><br><br><br><p style=\"right: 35px; top:100px; position:absolute;\">Display/Hide Stats: <input type= \"checkbox\" checked  onclick =\"cachedecache('LAND_STAT')\" ></p> <p id=\"LAND_STAT\" class =\"decache\">%s\n",st);
	

}

void save_SVG( struct landscape *fl,  char *file_out, AllParams *param,int R,char *st)
{
	char commande [1032];
	int n=(int)(strrchr(file_out,'.')-file_out);
	
//	draw_FL( 0, fl, file_out, 0,                  param->LAND_LOG,              1,                  -1,                    -1,                  0,                   1.0 ,R,  -1,  1,0,  -1, -1,0,0);
	draw_FL( 0, fl, file_out, param->LAND_THRESHOLD, param->LAND_LOG, !param->LAND_OPTCLEAN, param->LAND_DRAWFROM,param->LAND_DRAWTOEND, param->LAND_REFERENCE, param->LAND_SCALE_H/100.0,R, param->LAND_ONLY_MUT,1,param->LAND_CHAINS,1000,1000,param->LAND_COMPACT,param->LAND_FLAT,st);         

	sprintf(commande,"%s -f pdf -o %.*s.pdf %s",PROG_CONV,n,file_out,file_out);

	system(commande);
}

void save_Landscape ( struct landscape *fl,  char *Thefile)
{
	
	int g,i;
	int *geno=NULL;
	FILE *f;
	f=fopen(Thefile,"w");
	if (f==NULL)
		printf("pb...\n"),exit(1);
	
	for(i=0; i < fl->nlocus-1; i++)
	{
		fprintf(f,"%d ", (fl->alleles)[i] );
	}
	fprintf(f,"%d", (fl->alleles)[fl->nlocus-1] );
	fprintf(f,"\n");
	
	//int n=(int)(strrchr(file_fl,'.')-file_fl);
	
	for(g=0; g<fl->ngenotypes; g++ )
	{
		
		geno = int2genotype( *fl, g , geno );
		
		for(i=0; i < fl->nlocus-1; i++)
		{
			fprintf(f,"%d ", geno[i] );
		}
		fprintf(f,"%d", geno[fl->nlocus-1] );
		fprintf(f," %f\n",  fl->fitness[g] );
		
	}
	fclose(f);
	free(geno);
}

void outputparam(AllParams param,struct model_opt myoptions)
{
int i;
printf(" model %s<BR>",param.LAND_MODEL); // number associated to model see HTML order
printf(" model no %d<BR>",param.LAND_MODELNBR);
printf(" nb geno %d<BR>",param.LAND_NBLOCI);; //nbr of gentotype
printf(" same nbr %d<BR>",param.LAND_CHECK_SAMENBR);	// 0 if differents 1 if same
//printf(" nb interac %d",param.LAND_NK_PARAM); //if kaufman nbr of interactive loci
//printf(" rand %d",param.LAND_RAND_KAUF);	// 0 is not rand 1 if rand
printf(" log %d<BR>",param.LAND_LOG);
printf(" ref %d<BR>",param.LAND_REFERENCE); //nbr of left ref 
printf(" ratio fit %f<BR>",param.LAND_THRESHOLD); //ratio fitnees
printf(" clean %d<BR>",param.LAND_OPTCLEAN);// dsiplay neutrals
printf(" draw from %d<BR>",param.LAND_DRAWFROM); //display lnks from
printf(" missing %d<BR>",param.LAND_MISSING); //1 if 0 
printf(" w %f<BR>",param.LAND_SCALE_W); //width scale
printf(" h %f<BR>",param.LAND_SCALE_W); //heigth scale
printf(" nb fields %d<BR>",param.nbFields); // 
printf(" defaults %s<BR>", param.LAND_DEFAULT);
printf(" old =(%s)<BR>", param.LAND_PREVIOUS_FILE );
printf("Mult %f %f %f<BR>\n", myoptions.mu_s, myoptions.sigma_s, myoptions.DimRet);
printf("HoC %f<BR>\n",myoptions.sigma_hoc);
printf("EGG %f %f <BR>\n",myoptions.mu_eggbox,myoptions.sigma_eggbox);
printf("NK %d %d<BR>",myoptions.kaufman_K,myoptions.kaufman_rand);
printf("Ising %f %f %d <BR>\n",myoptions.mu_ising,myoptions.sigma_ising,myoptions.circular_ising);
printf("Opt %f %f %f %f <BR>\n",myoptions.mu_optimum,myoptions.mu_prod,myoptions.sigma_optimum, myoptions.sigma_prod);
printf("fixed= %f\n<BR>",myoptions.fix);
for (i=0;i<param.LAND_NBLOCI;i++)
	printf("%d,",param.ALLELS[i]);
printf("<BR>");
for (i=0;i<param.LAND_NBLOCI;i++)
	printf("%d,",param.LAND_MASK[i]);
printf("<BR>");	
}

void initparam(AllParams *param)
{
char *deft="0,1,1,-1,-1,0,1,1,1,0";
param->LAND_MODEL=NULL; //model name
param->LAND_NBLOCI=-1; //nbr of gentotype
param->LAND_MODELNBR=-1;
param->ALLELS=NULL; //how many allels in each genotype
param->LAND_CHECK_SAMENBR=0;	// 0 if differents 1 if same, if different the string is not present in the cGI field so initialise to 0
param->LAND_LOG=0;
param->LAND_REFERENCE=-1; //nbr of left ref 
param->LAND_THRESHOLD=1; //ratio fitness
param->LAND_OPTCLEAN=0;// dsiplay neutrals
param->LAND_DRAWFROM=-1; //display lnks from
param->LAND_DRAWTOEND=-1; //display lnks from
param->LAND_MISSING=0; //1 if misssing fitnesses are 0.0
param->LAND_SCALE_W=80; //width scale
param->LAND_SCALE_H=100; //heigth scale
param->nbFields=0; // 
param->LAND_DEFAULT=deft;
param->LAND_RMFWEIGHT=-1;
param->LAND_PREVIOUS_FILE =NULL;
param->LAND_MULTSAME=0;
param->LAND_CHAINS=0;
param->LAND_VALUE_MULTSAME=0;
param->LAND_ONLY_MUT=-1;

param->LAND_COMPACT=0;
param->LAND_FLAT	=0;
}

void freeparam(AllParams *param)
{
free(param->ALLELS);
free(param->LAND_MASK);
free(param);
}



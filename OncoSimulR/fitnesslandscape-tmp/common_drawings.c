

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>


#include "landscape.h"
#include "genotypes.h"
#include "drawings.h"

#define SCALE 0.1
#define COMPACT_SIZE 3
/*
	returns the nbr of mut from 00000 associated to the int x
*/
int dist_to_genoref( struct landscape *h, int *g, int *ref ){


	int l = h->nlocus-1;
	int nbmut=0;

	while( l >= 0 ){
		
		if(g[l] != ref[l])nbmut++;
		l--;
	}

	return nbmut;
}




/*
	Return a list of gentoypes with all genotypes that are neighbors and are >= divergent from 000...000
	fl is the fitness landscape ; g the current genotype
	li is memory for a list ; dist the distance of g to the 000 genotype 
*/
int RetrieveRightNeighbors( struct landscape *fl, int *genotype, int *geno_ref, struct list *li ){


	int x=0,          /* number of neighbor right genotypes */
	    d,            /* hamming distqnce to ref */
	    l,            /* locus */
	    a,            /* allele */
	    g2,           /* neighbor */
	    tmp;          /* used to built neighbors (g2) */
	

	d = dist_to_genoref( fl,  genotype , geno_ref );      /* hamming distance to ref */

	for( l=0; l < fl->nlocus; l++ ){
	
		tmp = genotype[l];
		
		for(a=0; a< fl->alleles[l]; a++){
		
			if( a != tmp ){ 
			
				genotype[l] = a;     /* a neihboring genotype */

 				if( dist_to_genoref( fl,  genotype , geno_ref ) > d || ( dist_to_genoref( fl,  genotype , geno_ref ) == d && a>tmp ) ){
										
					if( li->n < x+1 )
						resize_list( li->n+1, li );

					g2 = genotype2int( *fl, genotype);
					li->genotypes[x] = g2;
					
					x++;
				}
				
			}
		}
		genotype[l] = tmp;
	}

	return x;
}



/*
	bunch of svg commands rewritten to look the C way
*/
void SVG_header(FILE *f,int w, int h,char *s)
{
	fprintf(f,"<?xml version=\"1.0\" standalone=\"no\"?>\n");
	fprintf(f,"<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"\n");
	fprintf(f,"\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n");
	fprintf(f,"<svg xmlns=\"http://www.w3.org/2000/svg\"");
	fprintf(f," width=\"%d\" height=\"%d\">\n",w,h);
	fprintf(f,"<g>\n");
	fprintf(f,"<title>%s</title>\n",s);
}

void SVG_draw_line(FILE *f,int x_from,int y_from, int x_to,int y_to,char *color,int size)
{
	fprintf(f,"<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\"  style=\"stroke: %s; stroke-width:%f%%;\" />\n", x_from, y_from, x_to, y_to, color, size*SCALE);
}


void SVG_draw_dashline(FILE *f,int x1,int y1, int x2,int y2,char *color)
{
	fprintf(f,"<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\"  style=\"stroke: %s; stroke-width:%f%%; stroke-dasharray: 9, 5;\" />\n",x1,y1,x2,y2,color, SCALE);
}



void SVG_draw_text(FILE *f,int x1,int y1,char *s,int R)
{

	fprintf(f,"<text id=\"leg\" x=\"%d\" y=\"%d\" style=\"font-size:%dpt;\">%s</text>\n",x1,y1, 2*R ,s);
}

void SVG_draw_color_text(FILE *f,int x1,int y1,char *s,int R,char *color)
{

	fprintf(f,"<text id=\"leg\" x=\"%d\" y=\"%d\" style=\"font-size:%dpt; fill:%s;\">%s</text>\n",x1,y1, 2*R ,color,s);
}

void SVG_draw_float(FILE *f,int x1,int y1,float v,int size, char *color)
{
	fprintf(f,"<text x=\"%d\" y=\"%d\"  style=\"font-size:%dpt;fill:%s;\"  >%.3g</text>\n",x1,y1,size,color, v);		
}

void SVG_draw_int(FILE *f,int x1,int y1,int v,int size)
{
	fprintf(f,"<text x=\"%d\" y=\"%d\"  style=\"font-size:%dpt;\" >%d</text>\n",x1,y1,size,v);		
}


/*
	Specific SVG drawings
	draws a line of circles 1st one center is x,y
	g is the genotype id
*/
char *cold2warm( float frac, char *color ){

	unsigned char Red=0, Green=0, Blue=0;

	int c = (int) floor( frac*511 );

	
	if(color == NULL){
		color = (char *)calloc(20, sizeof(char));
		if(color == NULL)fprintf(stderr, "cold2warm: cannot allocate color, bye\n"), exit(3);
	}

	c=512+c;
	
	if(frac==1)
		c=1023;
		
	if(frac==0)
		c=512;
	
	if(c<256){
		Red=0;
		Blue=c;
		Green=255;
	}
	else if( c<512){
		Red=0;
		Blue=255;
		Green=511-c;
	}
	else if( c<768 ){
		Red=c-512;
		Blue=255;
		Green=100;
	
	}else{
		Red=255;
		Blue=1023-c;
		Green=100;
	}


	sprintf(color, "#%02X%02X%02X", Red,Green,Blue);

	return color;
}

/*
#define NBCOLORS 6
char *cold2warm( float frac ){
	

	char *LavaColors[NBCOLORS]={"#FFFF00","#FFCC00","#FF9900","FF6600","#FF3300","#993300"};
//	char *LavaColors[NBCOLORS]={"#CCFFFF","#CCCCFF","#CC99FF","#CC66FF","#CC33FF","#9900CC"};
//	char *LavaColors[NBCOLORS]={"#330066","#660066","#990066","CC0066","#FF0066","FF3300"};
	int c = (int) floor( frac*NBCOLORS);

//	printf("%f %d ",frac,c);

	if (c>NBCOLORS || c<0) printf("</svg> <p>PB in cold2warm %f %d</html>",frac,c),exit(1);
	
	return LavaColors[c];


}
 */
 
 
 void getcolor(float v,int *colr,float min,float max)
{


int v1;

colr[0]=0;
colr[1]=100;

v1=100*(max-v)/(max-min);

if (v1>50)
	colr[0]=(v1-50)*2; 
if (v1<50) //si v est proche de mid colr[1]=100 sinon tend vers 0...
	{
	colr[1]=v1*2;
	}
}

void draw_geno(int web,FILE *f, struct landscape *fl, int g, int **coord, char opt_log,int R,int from,int to,int compact,int flat)
{
	int i;

	int width=((fl->nlocus)* 2*R);
	int height=3*R	;
	char *locus_color[16]={"#c6c7c8","#1b171b","blue", "red",   "yellow", "green","purple", "lime","maroon", "navy",      /* color of the alleles */
					  "olive", "silver", "teal","fuchsia","aqua","gray"};
//0 is the lightest					  
//	char *depth_color[14]={"#AEEDDF","#7BE2E3","#50DBE8","#41D9ED","#31D4EC","#22CCEA","#1AB9E1","#12A8D4","#088EC7","#0876B6","#066FA2","#0A5E8D","#133670","#1B1740"};
//	char *depth_color[12]={"#1B1740","#133670","#0A5E8D","#066FA2","#0876B6","#088EC7","#12A8D4","A1AB9E1","#22CCEA","#50DBE8","#7BE2E3","#AEEDDF"};
//char *depth_color[12]={"#F8F8F8 ","#E8E8E8","#D8D8D8 ","#B8B8B8","#A8A8A8","#989898","#888888","#787878","#686868","#585858","#404040","#181818" }; //hight fit darkcolor
char *depth_color[12]={"#181818", "#404040", "#585858","#686868","#787878","#888888","#989898","#A8A8A8","#B8B8B8","#D8D8D8 ","#E8E8E8",  "#F8F8F8 " };// hight fit light color
	short c=0;

	char *colorf=NULL;       /* contain an RGB that encodes the warmth of the color -- from blue to red -- for genotype */
	int kol=0;
	int *genotype = NULL;    /* an int * version of genotype g */ 
	
	float frac=0;               /* how is fitness when compared to min and max */
	
	int nfitter;              /* number of fitter neighbor -- assess peak or well -- */
	int nworst;              /* number of fitter neighbor -- assess peak or well -- */

	
	
	if (compact)
		width=COMPACT_SIZE*R;
		
	genotype = int2genotype( *fl, g, genotype );
	
	/*
		Set how much relative fitness is high
	*/
	if(fl->maxf==fl->minf)
		kol=0;
	else
	{	
	if(opt_log)
		frac = (log( fl->fitness[g]) - log(fl->minf) ) / ( log(fl->maxf)-log(fl->minf) );
	else
		if (fl->maxf != fl->minf)
		frac = (fl->fitness[g] -fl->minf) / ( fl->maxf - fl->minf );
	kol=(int)(frac*11.0); 

	}	
		
	/*
		Set color warmth, based on frac
	*/
	colorf = cold2warm( frac, colorf );
	
//	colorf = cold2warm( frac);		
	/*
		Is a Peak or a Sink ?
	*/
	nfitter = CountFitterNeighbors( fl, g, 1, 0 );
	nworst  = fl->neighbors - CountFitterNeighbors( fl, g, 1, 1 );
	
	if( (nfitter==0 || nworst==0)  && CountDefinedNeighbor( fl, g) == fl->neighbors )
		c=1;


	/*
		Print the ellipse  if web is on put some javascript...
	*/
    if (web) 
    { 
    int x;
    int sizeoftriange=10;
    x=coord[0][g]+1 +(R*2*fl->nlocus)+2;
	// in case 1st genotype is a sink or a peak overwrite its color	
	if(c!=0)
		{
		
		if (nworst==0) //best fitness ever
			{
			if (!flat)
				{
				fprintf(f,"<g id=\"%d\" onclick=\"   displayID('%d','#e0141c',evt,1)\">\n",g,g);
				fprintf(f,"<rect id=\"%d\"   x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" rx=\"%d\" ry=\"%d\" style=\"fill:#e0141c;\" />\n", 
		          g,coord[0][g]-R, coord[1][g]-(height/2), width,height,3,3  );
		  		}
		    else
		       {  
		     
		    	fprintf(f,"<g id=\"%d\" onclick=\"   displayID('%d','#e0141c',evt,1)\">\n",g,g);
				fprintf(f,"<rect id=\"%d\"   x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" rx=\"%d\" ry=\"%d\" style=\"fill:%s;stroke:#e0141c;stroke-width:3;\" />\n", 
		          g,coord[0][g]-R, coord[1][g]-(height/2), width,height,3,3,depth_color[kol] );   
		        }   
		     }    
		else
			{
			if (!flat)
				{
				fprintf(f,"<g id=\"%d\" onclick=\" displayID('%d','#59ac2a',evt,1)\">\n",g,g);
				fprintf(f,"<rect id=\"%d\"   x=\"%d\" y=\"%d\"  width=\"%d\" height=\"%d\"  rx=\"%d\" ry=\"%d\" style=\"fill:#59ac2a;\" />\n", 
		           g,coord[0][g]-R, coord[1][g]-(height/2), width,height, 3,3 );
		     	}
		     else
		     	{
		     	fprintf(f,"<g id=\"%d\" onclick=\"   displayID('%d','#e0141c',evt,1)\">\n",g,g);
				fprintf(f,"<rect id=\"%d\"   x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" rx=\"%d\" ry=\"%d\" style=\"fill:%s;stroke:#59ac2a;stroke-width:3;\" />\n", 
		          g,coord[0][g]-R, coord[1][g]-(height/2), width,height,3,3,depth_color[kol] );   
		       }      
		    } 
		}           
	else
		{
		if (!compact)
			{
			 if (!flat)
		 		{
				fprintf(f,"<g id=\"%d\" onclick=\" displayID('%d','transparent',evt,0)\">\n",g,g);
				fprintf(f,"<rect id=\"%d\"   x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" rx=\"%d\" ry=\"%d\" style=\"fill:transparent;\" />\n", 
		            g,coord[0][g]-R, coord[1][g]-(height/2), width,height,3,3  );
		        }  	
		         else
		          
		     	{
//		     	fprintf(f,"<g id=\"%d\" onclick=\"displayIDdegrad('%d',evt)\">\n",g,g);
				fprintf(f,"<g id=\"%d\" onclick=\"displayID('%d','transparent',evt,0)\">\n",g,g);
		     	fprintf(f,"<rect id=\"%d\"   x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" rx=\"%d\" ry=\"%d\" style=\"fill:%s;\" />",
		     	g,coord[0][g]-R, coord[1][g]-(height/2), width,height,3,3 ,depth_color[kol]);
//				fprintf(f,"<rect id=\"%d\"   x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" rx=\"%d\" ry=\"%d\" style=\"fill:url(#MonDegrade%d);opacity:0.65\" />\n", 
//		         g,coord[0][g]-R, coord[1][g]-(height/2), width,height,3,3 ,g );
 		     	}  
		     }
		 else
		 	{
				fprintf(f,"<g id=\"%d\" onclick=\" displayID('%d','#cccccc',evt,0)\">\n",g,g);
				fprintf(f,"<rect id=\"%d\"   x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" rx=\"%d\" ry=\"%d\" style=\"fill:#cccccc;\" />\n", 
		            g,coord[0][g]-R, coord[1][g]-(height/2), width,height,3,3  );
		    }        
		}
	if (!compact)	
		{
/*		int rayon;
		if (opt_log)
		rayon=(R/2)+(2*R*(log(fl->fitness[g]))/(log(fl->maxf)-log(fl->minf)));
		else
		rayon=(R/2)+(2*R*(fl->fitness[g])/(fl->maxf-fl->minf));*/
//		if (!flat)
		for (i=0;i<fl->nlocus;i++)
			fprintf(f,"<circle cx=\"%d\" cy=\"%d\" r=\"%d\" style=\"  fill: %s\"   />\n", \
		            coord[0][g]+(i*2*R)+1,coord[1][g],2*R/3,  locus_color[(int)genotype[i]]);	
/*		else
		for (i=0;i<fl->nlocus;i++)
			fprintf(f,"<circle cx=\"%d\" cy=\"%d\" r=\"%d\" style=\"  fill: %s\"   />\n", \
		            coord[0][g]+(i*(rayon*2)+1),coord[1][g],rayon,  locus_color[(int)genotype[i]]);	*/
		 }           
	fprintf(f,"</g>");	           

	//if somes nodes where selected then add a triangle close to it	           
	if (from!=-1 && from==g)
		{
		if (compact)
			fprintf(f,"<path d=\"M %d %d l %d,-%d l 0,+%d l -%d,-%d \" stroke = \"#e0141c\" stroke-width = \"1\" fill = \"#e0141c\"/>\n", //red
				   coord[0][g]+width,coord[1][g],sizeoftriange,sizeoftriange,2*sizeoftriange,sizeoftriange,sizeoftriange);
		else
			fprintf(f,"<path d=\"M %d %d l %d,-%d l 0,+%d l -%d,-%d \" stroke = \"#e0141c\" stroke-width = \"1\" fill = \"#e0141c\"/>\n",
				   x,coord[1][g],sizeoftriange,sizeoftriange,2*sizeoftriange,sizeoftriange,sizeoftriange);
		}	

	if (to!=-1 && to==g)
		{
		if (compact)
			fprintf(f,"<path d=\"M %d %d l %d,-%d l 0,+%d l -%d,-%d \" stroke = \"#59ac2a\" stroke-width = \"1\" fill = \"#59ac2a\"/>\n",
				   coord[0][g]+width,coord[1][g],sizeoftriange,sizeoftriange,2*sizeoftriange,sizeoftriange,sizeoftriange);
		else
			fprintf(f,"<path d=\"M %d %d l %d,-%d l 0,+%d l -%d,-%d \" stroke = \"#59ac2a\" stroke-width = \"1\" fill = \"#59ac2a\"/>\n",
				   x,coord[1][g],sizeoftriange,sizeoftriange,2*sizeoftriange,sizeoftriange,sizeoftriange);
		}
	} //end of web
	else
	{    
	int x;
    int sizeoftriange=10;
    x=coord[0][g]+1 +(R*2*fl->nlocus)+2;
	// in case 1st genotype is a sink or a peak overwrite its color	
	if(c!=0)
		{
	 
		if (nworst==0) //sink
			{
			if (!flat)
			{
			fprintf(f,"<g id=\"%d\" >\n",g);
			fprintf(f,"<rect id=\"%d\"   x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" rx=\"%d\" ry=\"%d\" style=\"fill:#e0141c;\" />\n", 
		          g,coord[0][g]-R, coord[1][g]-(height/2), width,height,3,3  );
		     }
		     else{	fprintf(f,"<g id=\"%d\" >\n",g);			fprintf(f,"<rect id=\"%d\"   x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" rx=\"%d\" ry=\"%d\" style=\"fill:%s;stroke:#e0141c;stroke-width:3;\" />\n", 
		          g,coord[0][g]-R, coord[1][g]-(height/2), width,height,3,3,depth_color[kol] );   
}      
		    }       
		else
			{
			if(!flat)
			{
			fprintf(f,"<g id=\"%d\" >\n",g);
			fprintf(f,"<rect id=\"%d\"   x=\"%d\" y=\"%d\"  width=\"%d\" height=\"%d\"  rx=\"%d\" ry=\"%d\" style=\"fill:#59ac2a;\" />\n", 
		           g,coord[0][g]-R, coord[1][g]-(height/2), width,height, 3,3 );
		    }   
		     else{fprintf(f,"<g id=\"%d\" >\n",g);			fprintf(f,"<rect id=\"%d\"   x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" rx=\"%d\" ry=\"%d\" style=\"fill:%s;stroke:#59ac2a;stroke-width:3;\" />\n", 
		          g,coord[0][g]-R, coord[1][g]-(height/2), width,height,3,3,depth_color[kol] ); }        
		     }      
		}           
	else
		{
		
		if (!compact)
			{
			if (!flat)
			{
			fprintf(f,"<g id=\"%d\" >\n",g);
//			fprintf(f,"<rect id=\"%d\"   x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" rx=\"%d\" ry=\"%d\" style=\"fill:;\" />\n", 
//		           g,coord[0][g]-R, coord[1][g]-(height/2), width,height,3,3  );
			fprintf(f,"<rect id=\"%d\"   x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" rx=\"%d\" ry=\"%d\" style=\"fill:transparent;opacity:0\" />\n", 
		           g,coord[0][g]-R, coord[1][g]-(height/2), width,height,3,3  );
		     }
		    else
		    
		    	{
		     	fprintf(f,"<g id=\"%d\" >\n",g);
//				fprintf(f,"<rect id=\"%d\"   x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" rx=\"%d\" ry=\"%d\" style=\"fill:url(#MonDegrade%d);opacity:0.5\" />\n", 
//		         g,coord[0][g]-R, coord[1][g]-(height/2), width,height,3,3 ,g );
fprintf(f,"<rect id=\"%d\"   x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" rx=\"%d\" ry=\"%d\" style=\"fill:%s;\" />\n", 
		          g,coord[0][g]-R, coord[1][g]-(height/2), width,height,3,3,depth_color[kol] ); 
		     	}  

		    
		           
		     }
		 else
		 {
		
		 fprintf(f,"<g id=\"%d\" >\n",g);
		 fprintf(f,"<rect id=\"%d\"   x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" rx=\"%d\" ry=\"%d\" style=\"fill:#cccccc;\" />\n", 
		           g,coord[0][g]-R, coord[1][g]-(height/2), width,height,3,3  );
		    }        

		} //end of c==0
		
	if (!compact)	
		for (i=0;i<fl->nlocus;i++)
			fprintf(f,"<circle cx=\"%d\" cy=\"%d\" r=\"%d\" style=\"fill: %s;\"   />\n", 
		          coord[0][g]+(i*2*R)+1,coord[1][g],2*R/3,  locus_color[(int)genotype[i]]);	

	fprintf(f,"</g>");	           
	//if somes nodes where selected then add a triangle close to it	           
	if (from!=-1 && from==g) //draw red
	{
	if (!compact)	
		fprintf(f,"<path d=\"M %d %d l %d,-%d l 0,+%d l -%d,-%d \" stroke = \"#e0141c\" stroke-width = \"1\" fill = \"#e0141c\"/>",
				   x,coord[1][g],sizeoftriange,sizeoftriange,2*sizeoftriange,sizeoftriange,sizeoftriange);
	else
		fprintf(f,"<path d=\"M %d %d l %d,-%d l 0,+%d l -%d,-%d \" stroke = \"#e0141c\" stroke-width = \"1\" fill = \"#e0141c\"/>",
				   coord[0][g]+width,coord[1][g],sizeoftriange,sizeoftriange,2*sizeoftriange,sizeoftriange,sizeoftriange);
				   
	}			   
	if (to!=-1 && to==g)
		{
		if (!compact)	
		fprintf(f,"<path d=\"M %d %d l %d,-%d l 0,+%d l -%d,-%d \" stroke = \"#59ac2a\" stroke-width = \"1\" fill = \"#59ac2a\"/>",
				   x,coord[1][g],sizeoftriange,sizeoftriange,2*sizeoftriange,sizeoftriange,sizeoftriange);
		else
		fprintf(f,"<path d=\"M %d %d l %d,-%d l 0,+%d l -%d,-%d \" stroke = \"#59ac2a\" stroke-width = \"1\" fill = \"#59ac2a\"/>",
				   coord[0][g]+width,coord[1][g],sizeoftriange,sizeoftriange,2*sizeoftriange,sizeoftriange,sizeoftriange);
				   
		}		   
		
	} // end of !web
	
		
	
	
	
	
	free( genotype );
//	free(colorf);
}



/*
	draw a link between genotype g1 and g2
	ls is line size
*/
void draw_link(FILE *f, struct landscape *fl, int *geno1, int *geno2, int *geno_ref,  int **coord, int geno_width, char opt_neutral ,int R, float Fitness_Increase,int flat){
	char *colorn[3]={"#59ac2a","#e0141c","orange"};
//	char *colorbw[3]={"#50DBE8","#133670","lightgrey"}; //light blue darkblue
//	char *colorbw[3]={"#404040","#e0141c","lightgrey"}; //light grey darkgrey
	char *colorbw[3]={"#59ac2a","#e0141c","orange"}; //RDU
	
	char **color;
	
	if (!flat)
		color=colorn;
	else
		color=colorbw;
	
//	char *color[3]={"black","green","orange"};     /* color for gain, loss, neutral */
	int ls[3] = {1, 1, 1};                         /* line size */
	int c,cf;                                         /* indice of color */

	int g1 = genotype2int( *fl, geno1 ),            /* genotypes 1 and 2 and their hamming distance to ref */
	    g2 = genotype2int( *fl, geno2 ),
	    d1,
	    d2;                               


	int x1, x2, y1, y2;                            /* coordinates */

	/*
		set distances
	*/
	d1 = dist_to_genoref(fl, geno1, geno_ref );
	d2 = dist_to_genoref(fl, geno2, geno_ref );

	/*
		if not correctly ordered, re-order them g2 should be at the right of g1
	*/
	if(d2<d1){
		int tmp, *ptmp;
		tmp=g1; g1=g2; g2=tmp;
		tmp=d1; d1=d2; d2=tmp;
		ptmp=geno1; geno1=geno2; geno2=ptmp; 
	}


	/*
		set colors
	*/
	
	
	cf=  compare_fitness( fl -> fitness[g1], fl -> fitness[g2], Fitness_Increase );

	c = (fl -> fitness[g2] > fl -> fitness[g1])?0:1; 

//	if ( opt_neutral ) c = 2;

	

		
	x1 = coord[0][g1]+geno_width;
	x2 = coord[0][g2]-2*R;
	y1 = coord[1][g1];
	y2 = coord[1][g2];

	/*
		if equal, handle the special case --more than 2 alleles per locus--
	*/
	if(d1==d2){
		x1 -= (int) geno_width*4/5.0;
		x2 += (int) geno_width*4/5.0;
		y1+=R;
		y2-=R;
	}
	
	
//	fprintf(f,"<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\"  style=\"stroke: %s; stroke-width: %f%%; stroke-dasharray: 1, 1;\" />\n", x1,y1,x2,y2 ,color[c], ls[c]*SCALE);
	if (cf==0 && opt_neutral)
		fprintf(f,"<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\"  style=\"stroke: %s; stroke-width: %d; \" />\n", x1,y1,x2,y2 , color[2],ls[c]);
	else
		if (cf!=0)
			fprintf(f,"<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\"  style=\"stroke: %s; stroke-width: %d; \" />\n", x1,y1,x2,y2 ,color[c], ls[c]);
		
}

#define SIGN(a) ( ( (a)<0 )?-1:1 )

// renvoie 1 si neutre
int isneutral(float f1, float f2,float threshold)
{
float ratio;
			if(f1 == 0 && f2 == 0 )
				ratio = 1;
			else if(f1 == 0 || f2 == 0 || SIGN(f1) != SIGN(f2))
				ratio = 1000000;
			else
				ratio = ( fabs(f2) > fabs(f1) ) ? f2/f1 : f1/f2;
				
if (ratio <=threshold)
return 1;
else
return 0;
}

/*
	Draw lines between points coord 
	
	(x,y) in **coord are the points at the begining of geno if from this have to add decallage to the 1st point
	line is dashline when values are under threshold, red when there is a fitness decrease,
	green if stable and black otherwise
	OnlyMutation affiche uniquement la mutation passée default -1
	the larger the line is, the bigger the fitness increase
*/

void draw_links( FILE *f, struct landscape *fl, int **coord, int geno_width, float threshold, int opt_log, int opt_clear, int *geno_ref,int R, int Only_Mutation ,int flat)
{
	
	
	int g,                 /* current genotype */
	    g2,                /* another genotype (i.e. here it is a right neighbor) */
	    i,j;
	    
	struct list li_rn;     /* list of right_neighbors */     
	int n_rn=0;            /* number of right neighbor with distance to ref >= d */
	
	int *genotype=NULL;    /* int * version of genotype */
	int *geno2=NULL;       /* int * version of another genotype */


	threshold = (threshold<1)?1:threshold;   /* ratio can never be smaller than 1 */

	li_rn = new_empty_list();
		
	
	/*
		For each genotypes
			only draw link to right neighbor (i.e. with a higher distance to ref)
	*/

	for (g=0;g<fl->ngenotypes;g++)
	{
	
		if( fl->fitness[g] == DEFAULT_FITNESS )
			continue;
		
		genotype = int2genotype( *fl, g, genotype);                               /* genotype in int */
		n_rn     = RetrieveRightNeighbors( fl, genotype, geno_ref, &li_rn  );    /* all neighbors with distance to ref >=d (n_rn is their #, li_rn their list) */
 		
		/*
			scan list of right neighbors 
		*/
		for (i=0;i<n_rn;i++)
		{
			int mutation = -1;

			g2 = li_rn.genotypes[i];
			geno2 = int2genotype( *fl, g2, geno2 );

			for( j=0;j<fl->nlocus;j++ )
				if( geno2[j] != genotype[j] )
					mutation = j+1;


			if( fl->fitness[g2] == DEFAULT_FITNESS )
				continue;
		
	

			if ((isneutral(fl->fitness[g],fl->fitness[g2],threshold)==0) && ( Only_Mutation == -1 || (Only_Mutation != -1 && Only_Mutation == mutation) ) )
				{
				draw_link(f, fl, genotype, geno2, geno_ref, coord, geno_width,  0 , R ,threshold,flat);
				}
			else
				if( opt_clear == 0 && ( Only_Mutation == -1 || (Only_Mutation != -1 && Only_Mutation == mutation) ) )
					{
					draw_link(f, fl,  genotype, geno2, geno_ref, coord, geno_width, 1 , R,threshold ,flat);
					
					}
				else	/*SOPHIE opt_clear == 0 && */
				if (( opt_clear == 0 )&& isneutral(fl->fitness[g],fl->fitness[g2],threshold))
					{
					draw_link(f, fl,  genotype, geno2, geno_ref, coord, geno_width, 1 , R,threshold ,flat);
					}
						
		}
	}
	
	free(genotype);
	free(geno2);
	free_list(&li_rn);
}


/*Draw chains i.e. only nodes from where only one better link is present
can be used with limits from , to or both*/
void draw_to_chains( FILE *f, struct landscape *fl, int **coord, int geno_width,float threshold, int to, int *geno_ref  ,int R,int flat)
{

	struct list li = new_empty_list();
	struct list li2 = new_empty_list();
	int i;
	int g2;
	int x,y;

	int *geno=NULL;        /* int * version of the genotype */
	int *geno2=NULL;       /* int * version of another genotype */
	
	geno = int2genotype( *fl, to, geno );
	
	x = RetrieveNeighbors( fl,to, NULL, &li, "w", threshold );




	for(i=0;i<x;i++){
	
		g2 = li.genotypes[i];
		geno2 = int2genotype( *fl, g2, geno2 );
		y= RetrieveNeighbors( fl, g2, NULL, &li2, "f", threshold );
		if (y==1)
			{
			
			if( fl->fitness[g2] == DEFAULT_FITNESS )
				continue;
				
			draw_link(f, fl,  geno2, geno, geno_ref, coord, geno_width, 0 , R,threshold,flat);
				
		 	draw_to_chains( f, fl, coord, geno_width, threshold, g2, geno_ref ,R,flat);	
		 	}
	}
	
	
	free(geno);
	free(geno2);
	
	free_list(&li);
	free_list(&li2);


}








/*Draw chains i.e. only nodes from where only one better link is present
can be used with limits from , to or both*/

void draw_chains( FILE *f, struct landscape *fl, int **coord, int geno_width,float threshold, int from, int to, int *geno_ref  ,int R,int flat)
{
	struct list li = new_empty_list();
	int i,x,depart=-1,fin=-1,depth=0,flag,j;
	int *geno=NULL;        /* int * version of the genotype */
	int *geno2=NULL;       /* int * version of another genotype */
	int *chemin;
	chemin=malloc(sizeof(int)*2*fl->ngenotypes);
	
	// scan all genotype for nodes having only one best neighbour 
	//if from or to specified , look if they are found in such nodes mark debut and fin and store genotypes linked in chemin
	//therefore chemin is [a,b,b,c] for  a--b--c
	//links will be drawn after
	
		for (i=0;i<fl->ngenotypes;i++)
			{
			li=new_empty_list();
			
			x = RetrieveNeighbors( fl, i, NULL, &li, "f", threshold);
			geno = int2genotype( *fl, i, geno );
		
		
			if (x==1)
				{
				if (i==from)depart=depth;
				if (li.genotypes[0]==to)fin=depth+1;
				geno2=int2genotype( *fl, li.genotypes[0], geno2 );
				chemin[depth++]=i;
				chemin[depth++]=li.genotypes[0];
				if (from==-1 && to==-1)
					draw_link(f, fl,  geno, geno2, geno_ref, coord, geno_width, 0 , R,threshold,flat);
	
		
				}
			free_list(&li);	
			}
	
	if (to==-1 && (from!=-1 && depart!=-1)) //from is the begining of a chain draw no end specified
		{
		i=depart;
		while (1)
			{
			flag=0;
			geno = int2genotype( *fl, chemin[i], geno );
			geno2 = int2genotype( *fl, chemin[i+1], geno2 );

			draw_link(f, fl,  geno, geno2, geno_ref, coord, geno_width, 0 , R,threshold,flat);
			if (chemin[i+1]==fl->ngenotypes-1) break;
			for (j=i+2;j<depth;j++) //looking in chemin if there is a path from chemin[i+1]
				if (chemin[j]==chemin[i+1]) {flag=1;i=j;break;}
			if (flag==0)break;	// no path found so end of road
		
			}
		}		
	else
		if (from==-1 && (to!=-1 && fin!=-1)) //end is the end and no from
			{
			draw_to_chains( f, fl, coord, geno_width,threshold, to, geno_ref  , R,flat);
			
			}		
		else
			if (depart!=-1 && from!=-1 && to!=-1 && fin!=-1) //end and from specified and both find		
				{
				i=depart;
				while  (1)
					{
					flag=0;
					geno = int2genotype( *fl, chemin[i], geno );
					geno2 = int2genotype( *fl, chemin[i+1], geno2 );
					
					draw_link(f, fl,  geno, geno2, geno_ref, coord, geno_width, 0 , R,threshold,flat);
					
					if (chemin[i+1]== to) break; 
					for (j=i+2;j<depth;j++)
						if (chemin[j]==chemin[i+1]) {flag=1;i=j;break;}
					if (flag==0)break;	
					}
			
				}	
			
	free(chemin);
	
	free(geno);
	
	free(geno2);
	


	


}


/*
	Draw all links from a starting genotype
	
	(x,y) in **coord are the points at the begining of geno if from this have to add decallage to the 1st point
	line is dashline when values are under threshold, red when there is a fitness decrease,
	green if stable and black otherwise
	
	the larger the line is, the bigger the fitness increase
*/
	
void link_fitter_neighbor( FILE *f, struct landscape *fl, int **coord, int geno_width,float threshold, int g, int *geno_ref  ,int R,int opt_clear,char *visited_genotypes,int flat){

	struct list li = new_empty_list();
	struct list li2 = new_empty_list();
	int i;
	int g2;
	int x,y;

	int *geno=NULL;        /* int * version of the genotype */
	int *geno2=NULL;       /* int * version of another genotype */
	
	geno = int2genotype( *fl, g, geno );
	
	x = RetrieveNeighbors( fl, g, NULL, &li, "fn", threshold );

	for(i=0;i<x;i++){
		
		/* connect */
		g2 = li.genotypes[i];
		geno2 = int2genotype( *fl, g2, geno2 );
		
		if( fl->fitness[g2] == DEFAULT_FITNESS )
			continue;

		if( opt_clear == 0 )
			draw_link(f, fl,  geno, geno2, geno_ref, coord, geno_width, 1 , R,threshold,flat);
		else
			draw_link(f, fl,  geno, geno2, geno_ref, coord, geno_width, 0 , R,threshold,flat);
			
		/*draw_link( f, fl, geno, geno2, geno_ref, coord, geno_width, isneutral , opt_clear ,R);*/
	}
	
	if(opt_clear)
		y = RetrieveNeighbors( fl, g, visited_genotypes, &li2, "f", threshold );
	else
		y = RetrieveNeighbors( fl, g, visited_genotypes, &li2, "fn", threshold );

	for(i=0;i<y;i++){
	
		g2 = li2.genotypes[i];
		
		if( fl->fitness[g2] == DEFAULT_FITNESS )
			continue;
		
		/* call the function with g2 */
		link_fitter_neighbor( f, fl, coord, geno_width, threshold, g2, geno_ref ,R, opt_clear,visited_genotypes,flat);
		
	}
	
	free(geno);
	free(geno2);
	
	free_list(&li);
	free_list(&li2);

}



	
void link_cripple_neighbor( FILE *f, struct landscape *fl, int **coord, int geno_width,float threshold, int g, int *geno_ref  ,int R,int opt_clear,char *visited_genotypes,int flat){

	struct list li = new_empty_list();
	struct list li2 = new_empty_list();
	int i;
	int g2;
	int x,y;

	int *geno=NULL;        /* int * version of the genotype */
	int *geno2=NULL;       /* int * version of another genotype */
	
	geno = int2genotype( *fl, g, geno );
	
	x = RetrieveNeighbors( fl, g, NULL, &li, "wn", threshold );

	for(i=0;i<x;i++){
		
		/* connect */
		g2 = li.genotypes[i];
		geno2 = int2genotype( *fl, g2, geno2 );
		
		if( fl->fitness[g2] == DEFAULT_FITNESS )
			continue;

		if( opt_clear == 0 )
			draw_link(f, fl,  geno, geno2, geno_ref, coord, geno_width, 1 , R,threshold,flat);
		else
			draw_link(f, fl,  geno, geno2, geno_ref, coord, geno_width, 0 , R,threshold,flat);
			
		/*draw_link( f, fl, geno, geno2, geno_ref, coord, geno_width, isneutral , opt_clear ,R);*/
	}
	
	if(opt_clear)
		y = RetrieveNeighbors( fl, g, visited_genotypes, &li2, "w", threshold );
	else
		y = RetrieveNeighbors( fl, g, visited_genotypes, &li2, "wn", threshold );


	for(i=0;i<y;i++){
	
		g2 = li2.genotypes[i];
		
		if( fl->fitness[g2] == DEFAULT_FITNESS )
			continue;
		
		/* call the function with g2 */
		link_cripple_neighbor( f, fl, coord, geno_width, threshold, g2, geno_ref ,R, opt_clear,visited_genotypes,flat);
		
	}
	
	free(geno);
	free(geno2);
	
	free_list(&li);
	free_list(&li2);

}

		//	link_from_to_neighbor(  FILE *f, struct landscape *fl, int **coord, intwidth_geno,float threshold,int opt_from,int opt_to,int * geno_ref, R ,0 ,chemin, opt_clear,visited_links,nodes);
/*
static void link_from_to_neighbor( FILE *f, struct landscape *fl, int **coord, int geno_width,float threshold, int from,  int to, int *geno_ref  ,int R,int previous,int opt_clear,char *visited_links, char *nodes) {
	struct list li = new_empty_list();

	int i,n,g2,m;

	int *geno=NULL;       
	int *geno2=NULL;       

	if (from==to )    //yes found
		{	
			geno=int2genotype( *fl,previous,geno);
			geno2=int2genotype( *fl,to,geno2);
			visited_links[previous*fl->ngenotypes+to]=visited_links[to*fl->ngenotypes+previous]=1;
			nodes[previous]=1; //this node is part of a road
			if( opt_clear == 0 )
					draw_link(f, fl,  geno, geno2, geno_ref, coord, geno_width, 1 , R,threshold);
			else
					draw_link(f, fl,  geno, geno2, geno_ref, coord, geno_width, 0 , R,threshold);		
				
			return;
			}
	else
	if(	nodes[previous]==1	)
		{
			geno=int2genotype( *fl,previous,geno);
			geno2=int2genotype( *fl,from,geno2);
			visited_links[previous*fl->ngenotypes+from]=visited_links[from*fl->ngenotypes+previous]=1;
			nodes[previous]=1; //this node is part of a road
			if( opt_clear == 0 )
					draw_link(f, fl,  geno, geno2, geno_ref, coord, geno_width, 1 , R,threshold);
			else
					draw_link(f, fl,  geno, geno2, geno_ref, coord, geno_width, 0 , R,threshold);		
				
			return;
	
		}
		



	if(opt_clear)
		n = RetrieveNeighbors( fl, from, NULL, &li, "f", threshold );
	else
	    n= RetrieveNeighbors( fl, from, NULL, &li, "fn", threshold );


 	for( i=0;i<n;i++) {
		g2 = li.genotypes[i];
		
		link_from_to_neighbor(f,fl,coord,geno_width,threshold,g2,to,geno_ref,R, from , opt_clear,visited_links,nodes);
	}
 
	free_list(&li);
}
*/
#define HH 60
#define NTICS 4

static int link_from_to_neighbor( FILE *f, struct landscape *fl, int **coord, int geno_width,float threshold, int from,  int to, int *geno_ref  ,int R,int depth,int *chemin,int opt_clear,char *visited_genotypes,int flat) {
	struct list li = new_empty_list();

	int i,n,g2,c=0;

	int *geno=NULL;       
	int *geno2=NULL;       
	

	if(depth > fl->ngenotypes){printf("aaaargh why so much hate\n");exit (1);	}
	
	if(visited_genotypes[from]==1 || visited_genotypes[from]==-1)
		return -1000000;
	
	if( visited_genotypes[from]==0 )
		visited_genotypes[from] = 1;

	chemin[depth]=from;

	if (visited_genotypes[from]==2)
		{
		for( i=depth; i>0 ;i-- ) {
			
			geno=int2genotype( *fl,chemin[i-1],geno);
			geno2=int2genotype( *fl,chemin[i],geno2);
			
			
			//if( fl->fitness[chemin[i]] == DEFAULT_FITNESS )
			//	continue;

			if( opt_clear == 0 )
					draw_link(f, fl,  geno, geno2, geno_ref, coord, geno_width, 1 , R,threshold,flat);
			else
					draw_link(f, fl,  geno, geno2, geno_ref, coord, geno_width, 0 , R,threshold,flat);		
			
			if( visited_genotypes[chemin[i]]+visited_genotypes[chemin[i-1]] == 4 )
				break;

			visited_genotypes[chemin[i-1]]=2;

		 	
			}
		return 1;
	}


	if(opt_clear)
		n = RetrieveNeighbors( fl, from, NULL, &li, "f", threshold );
	else
	    n= RetrieveNeighbors( fl, from, NULL, &li, "fn", threshold );


 	for( i=0;i<n;i++) {
 	
 	
		g2 = li.genotypes[i];
		
		if (visited_genotypes[g2]!=-1)
			{
			c += link_from_to_neighbor(f,fl,coord,geno_width,threshold,g2,to,geno_ref,R,depth+1,chemin, opt_clear,visited_genotypes,flat);
			}
	}

	if (c==0)
 		visited_genotypes[from]=-1;
	
	
	free_list(&li);

	return c;
}
 
int maxcomb(struct landscape land, int *l)
{
	int maxi=0;
	int i;
	int c,*g;
	g=NULL;
	for (i=0;i<land.nlocus+1;i++)
		l[i]=0;

	for (i=0;i<land.ngenotypes;i++)
	{
	
		g=int2genotype( land, i, g);
	
		c=int2nbmut(land,g); 
	
		l[c]++;
	}

	for (i=0;i<land.nlocus+1;i++)
		if (maxi<l[i]) maxi=l[i];	
	
	return(maxi);		
}

int int2nbmut( struct landscape h, int *g ){


        int l=h.nlocus-1;
        int nbmut=0;



        while( l >= 0 ){

                if(g[l]>0)nbmut++;
                l--;
        }

        return nbmut;
}

/*enleves les espaces blancs ou tab en trop*/
int trim(char *l)
	{
	int n=strlen(l),i,j;
	char *newl;
	for (i=n-1;i>0;)
		{
		if (l[i]!=' ' && l[i]!='\t')
			break;
		i--;
		}	
	l[i+1]='\0'		;
	newl=malloc(sizeof(char)*(n+1));
	for (i=1,j=0;i<n;i++)
		{
		if ((l[i]==' '|| l[i]=='\t') && (l[i]==l[i-1]))
		   continue;
		else
		newl[j++]=l[i-1];   
		}
	
	
	if (l[i-1]!=' ' && l[i]!='\t')
	   	newl[j++]=l[i-1];   
	 newl[j]='\0';
	 
	 strcpy(l,newl);
	 free(newl);
	
	 j=0;n=0;
	while (l[j]!='\0')
		{
			if (l[j]==' ' || l[j]=='\t') n++;
			j++;
		}	

	return(n+1);

	}   	


    


void  parseHTML(char *stats)
{
char *p=stats,*m;
int i;

m=strstr(p,"<BR>");
while (m!=NULL)
{
for(i=0;i<3;i++)
	m[i]=' '; 
m[i]='?';
m=strstr(m+1,"<BR>");
}
m=strstr(p,"&nbsp;");
while (m!=NULL)
{
for(i=0;i<6;i++)
	m[i]=' '; 

m=strstr(m+1,"&nbsp;");
}
trim(stats);

}

void html_to_SVG(FILE *f,char *stats,int x,int y, int w)
{
int l,i=0,j=0;
parseHTML(stats);
l=strlen(stats);
while(1)
	{
	int c=0;
	while (i<l && (*(stats+i)!='?') && (*(stats+i)!='\0'))
	  	{i++;c++;}
	if (i>=l || (*(stats+i)=='\0'))  break;
  
	fprintf(f,"<text  x=\"%d\" y=\"%d\" style=\"font-size:12pt;\">%.*s</text>\n",x,y,c,stats+j);
	y+=15;
	
	i++;
	j=i;
	}
}
void get_coords(int *width,int *height,	int *coord[2], struct landscape *fl, int *geno_ref,int opt_log,int margin_w,int margin_h,int width_geno,int R,int compact)
{
	int g;                             /* the genotypes */
	int *genotype=NULL;                /* the int * version of g */
	int d;                             /* hamming distance to GenoRef */
	int *stack;                        /* how many genotypes at dist i have been stored at height j. i \in [0, nlocus+1], j is [0,49]*/
	float frac=0;                        /* where fitness fall between max and min --eventually with log-- */
	int new_width=0,new_height=0;
		


	stack = (int *)calloc( (size_t) (HH+1)*( fl->nlocus+1 ), (size_t)sizeof(int) );
//if( !stack )fprintf(ferr, "draw_FL: cannot allocate stack, bye\n"), exit(3);
	
	for (g=0 ; g<fl->ngenotypes ; g++)
	{
		genotype = int2genotype( *fl, g, genotype);
		d = dist_to_genoref ( fl, genotype , geno_ref); 
		
		if (fl->fitness[g] != DEFAULT_FITNESS){
			
			int bin=0;   /* bin number is between 0 and HH-1 */
			//float v;
			
			if(opt_log) 
				{
				if (fl->minf!=0)
				frac = (log( fl->fitness[g]) - log(fl->minf) ) / ( log(fl->maxf)-log(fl->minf) );
				else
				frac = (log( fl->fitness[g]) - log(fl->minf) ) /  log(fl->maxf);
				}
			else
				if (fl->maxf != fl->minf)
				frac = (fl->fitness[g] -fl->minf) / ( fl->maxf - fl->minf );
				
			
			/*
				check whether a genotype is already here, in that case, one must add some values to x
			*/

			bin = (int) lround( frac*HH );

			if(bin == HH)
				bin--;


			if(compact)	
			coord[0][g] = margin_w + (d * width_geno * 3) ;
			else	
			coord[0][g] = margin_w + (d * width_geno * 3) + stack[ d*HH + bin ] * 3*R/2;
			if(compact)	
			coord[1][g] = margin_h + (*height)*(1 - frac) ;
			else
			coord[1][g] = margin_h + (*height)*(1 - frac) -R + stack[ d*HH + bin ] * 3*R/2;     //+Stack              /* from 0,0 bottom left to 0,0 on top, left */
		
			
			
			if(new_width <coord[0][g])
				new_width =coord[0][g];
//			stack[ d*HH + (int) lround( frac*HH ) ] ++;
			stack[ d*HH + bin ] ++;
			if(new_height <coord[1][g])
				new_height =coord[1][g];
				
			
		}
	}
	
//Now resize if needed because when too much geno are in same position it goes outside of graph...


if (new_width >*width)
	*width=new_width;

	
if (new_height >*height)
	*height= new_height;

free (stack);
}

void get_coords_flat(int *width,int *height,	int *coord[2], struct landscape *fl,int margin_w, int margin_h, int width_geno,int R,int *geno_ref)
{
	int d,*g;

	int x,y,i;
	int  maxi_comb;
	int *l;
	int sp_h=4;
	int dec;

	int height_geno = 2*R;

	g=NULL;
	
	l=malloc(sizeof(int)*(fl->nlocus+1));
	
	

	maxi_comb=maxcomb(*fl,l);// how many geno max
	

	*height=(maxi_comb*sp_h*height_geno)+(2*R)+(2* margin_h);
	
	
	for (i=0;i<2;i++)
		coord[i]=malloc (sizeof(int)*fl->ngenotypes);

	for (i=0;i<fl->nlocus+1;i++)
	{

		dec=0;
		dec=(ceil(maxi_comb-l[i])/2.0)*sp_h*height_geno;
		l[i]= margin_h+R+dec;

	}	

	for (i=0;i<fl->ngenotypes;i++)
		{
		
		
		
		
		g=int2genotype( *fl, i, g);
//		c=int2nbmut(*fl,g); 
		d = dist_to_genoref ( fl, g , geno_ref); 
		
		x= margin_w+R+(d*width_geno*3);
		
		
		y=l[d];
		coord[0][i]=x;
		coord[1][i]=y;
	
		l[d]+=(sp_h*height_geno);
		}
		
		
free(l);	

	
}

/* draw_FL list of arguments
web=1 if draw web page
fl all the structure
fileout file ou stdout
threshold= value for determining neutrality
opt_log =1 if log scale
opt_clear=1 keep neutral links
opt_ref=begining of landscape, all links issued from this node and better are shown
opt_to=end of landscape, all links issued from this node and lower are shown
genoref= most right genotype
rescale= if changing scale
R =for scale also
Only_mutation=-1 for me
with_legen=1 if have to draw a legend
*/

void draw_FL( int web, struct landscape *fl,  char *file_out,  float threshold,  int opt_log,  int opt_clear, int opt_from, int opt_to, int genoref, float rescale_height ,int R, int Only_Mutation, int with_legend,int opt_chains,int wwii,int hhee,int compact,int flat,char *stats)
{
	FILE *f,*ferr;

	int i;                             /* dummy counters */

	int x=0,y=0;
	char *colorn[3]={"#59ac2a","#e0141c","orange" };   
	
	char legend[3][32]={"gain","loss","neutral"};
//	char *colorbw[3]={"#50DBE8","#133670","lightgrey"}; 
	char *colorbw[3]={"#404040","#e0141c","lightgrey"}; //light grey darkgrey

	int width,                         /* dimensions  */
	    height,
	    bout=0	;       	/*size of legend*/
	//int new_width=0,new_height=0,widthg,heighthg;
	int margin_w = 12*R;                  /* magins on height and width */
	int margin_h = 5*R;
	int add_leg; 						/*add room for legend*/
	int width_geno  = R*2*fl->nlocus;   /* dimension of a genotype */
	int height_geno = 2*R;
	int wstats=180;
//	int *stack;                        /* how many genotypes at dist i have been stored at height j. i \in [0, nlocus+1], j is [0,49]*/

	int *coord[2];                     /* coords of all genotypes in svg sheet */

	int g;                             /* the genotypes */
	int *genotype=NULL;                /* the int * version of g */
	int *geno_ref=NULL;                /* the int * version of genoref */

//	float frac;                        /* where fitness fall between max and min --eventually with log-- */

	float tics[NTICS];                 /* fitness values for the tics */
	int *chemin;
	char **color;
	int add_width=0,add_height=0;
	if (!flat)
		color=colorn;
	else
		color=colorbw;
	if (web)
	{f=stdout;ferr=stdout;}
	else
	{	
	f = fopen( file_out, "w" );
	ferr=stderr;
	if(!f)fprintf(stderr, "problem with output file, bye\n"), exit(4);
	}
	
	sprintf(legend[2],"%s (%.1f)",legend[2],threshold);

		

	width     = width_geno*((fl->nlocus+1)*3-2);
	height    = HH*height_geno*rescale_height;
//	fprintf(f,"%d %d<BR>",width ,height );exit(1);
	
		/*
		Compute and store all coordinates of all genotypes 
	*/
	
//		if (compact) width_geno  = R*2;

	coord[0] = (int *)malloc( (size_t) sizeof(int)*fl->ngenotypes );
	coord[1] = (int *)malloc( (size_t) sizeof(int)*fl->ngenotypes );
	if( !coord[0] || !coord[1] )fprintf(stderr, "draw_FL: cannot allocate coord, bye\n"), exit(3);

	for(g=0; g<fl->ngenotypes ; g++)
		coord[0][g]= coord[1][g] = 0;


	geno_ref = int2genotype( *fl, genoref, geno_ref );       /* the int* version of the reference --make it easy for hamming distance */
			
	if (!flat)
		get_coords(&width, &height,	coord, fl, geno_ref,opt_log,margin_w,margin_h,width_geno,R,compact);
	else
	    get_coords_flat(&width,&height,coord, fl,margin_w,margin_h,width_geno, R, geno_ref);
	
	if (with_legend)
		{
		add_leg=100;
		}
	else
		add_leg=0;
	// add some room for very small landscape...	
	if ((width + 2*margin_w)	<600)
		add_width=600-(width + 2*margin_w);
	if ((height + add_leg+ 2*margin_h<700))
		add_height=700-height + add_leg+ 2*margin_h;
	if (web)
		{

		fprintf(f,"	<div id=\"svg_content\">\n");
		fprintf(f,"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%d\" height=\"%d\">\n",(add_width+width + 2*margin_w), height + add_leg+ 2*margin_h);

		fprintf(f,"<g>\n");
		if (flat) //then define some degrades
			{
			fprintf(f,"  <defs>\n");
			int colr[2];
			char pc='%';
			for (g=fl->ngenotypes-1;g>=0;g--)
				{
				
			    fprintf(f,"<linearGradient id=\"MonDegrade%d\">\n",g);
			    if (opt_log)
			    getcolor(log(fl->fitness[g]),colr,log(fl->minf),log(fl->maxf));
			    else
    			getcolor(fl->fitness[g],colr,fl->minf,fl->maxf);
  				fprintf(f,"      <stop offset=\"%d%c\" stop-color=\"#e0141c\" />\n",colr[0],pc);
   	 			fprintf(f,"      <stop offset=\"%d%c\" stop-color=\"#59ac2a\" />\n",colr[1],pc);
    			fprintf(f,"  </linearGradient>\n");
    			}

			fprintf(f,"  </defs>\n");
			}
		}
	else
		{
		SVG_header(f,width+add_width+ 2*margin_w +wstats, add_height+height+ add_leg + 2*margin_h+bout, file_out);
		if (flat) //then define some degrades
			{
			fprintf(f,"  <defs>\n");
			int colr[2];
			char pc='%';
			for (g=fl->ngenotypes-1;g>=0;g--)
				{
				
			    fprintf(f,"<linearGradient id=\"MonDegrade%d\">\n",g);
			    if (opt_log)
			    getcolor(log(fl->fitness[g]),colr,log(fl->minf),log(fl->maxf));
			    else
    			getcolor(fl->fitness[g],colr,fl->minf,fl->maxf);
  				fprintf(f,"      <stop offset=\"%d%c\" stop-color=\"#e0141c\" />\n",colr[0],pc);
   	 			fprintf(f,"      <stop offset=\"%d%c\" stop-color=\"#59ac2a\" />\n",colr[1],pc);
    			fprintf(f,"  </linearGradient>\n");
    			}

			fprintf(f,"  </defs>\n");
			}
		}
	/*
		Set tics values -- NTICS tics for now --
	*/
	if (!flat) //draw axes if not flat
		{
		tics[0] = fl->minf;
		tics[NTICS-1] = fl->maxf;

		for(i=0 ; i<NTICS; i++)
			if(opt_log)
				{
				if (fl->minf!=0)
				tics[ i ] = exp(  log(fl->minf) + i*(log(fl->maxf)-log(fl->minf))/(NTICS-1.0)  );
				else
				tics[ i ] = exp(   i*(log(fl->maxf))/(NTICS-1.0)  );
				}
			else
				tics[ i ] = fl->minf + i*(fl->maxf-fl->minf)/(NTICS-1.0);
	
	
		/*
			Y-axis with fitness values
		*/
		SVG_draw_line(f, margin_w,       0, margin_w,       height+2*margin_h, "black", 1);
		SVG_draw_line(f, margin_w+width, 0, margin_w+width, height+2*margin_h, "black", 1);
	
		for(i=0; i<NTICS;i++){

			SVG_draw_line(f,  margin_w+width+10,  margin_h + (int)(( 1.0 - i/(NTICS-1.00) )*height),  margin_w+width,   margin_h + (int)(( 1.0 - i/(NTICS-1.00) )*height), "black", 1);
			SVG_draw_float(f, margin_w+width+10,  margin_h + (int)(( 1.0 - i/(NTICS-1.00) )*height),  tics[i],  2*R  ,"black");
			SVG_draw_line(f,  margin_w-10,     margin_h + (int)(( 1.0 - i/(NTICS-1.00) )*height),  margin_w,   margin_h + (int)(( 1.0 - i/(NTICS-1.00) )*height), "black", 1);
			SVG_draw_float(f, 1,  margin_h + (int)(( 1.0 - i/(NTICS-1.00) )*height),  tics[i],  2*R  ,"black");

			}
		}



	/*
		Connect the ellipses
	*/
if (compact) width_geno  = R*COMPACT_SIZE;

	chemin=malloc(sizeof(int)*fl->ngenotypes); //to store paths 
	if (opt_chains==1)
		draw_chains( f, fl, coord, width_geno,threshold, opt_from, opt_to, geno_ref  , R,flat);
	else	
		if(opt_from==-1 && opt_to== -1) // draw the whole landscape
			draw_links(f, fl, coord, width_geno, threshold, opt_log, opt_clear, geno_ref, R, Only_Mutation,flat);
		else

			{
			char *visited_genotypes;
			visited_genotypes = (char *)calloc( (size_t) fl->ngenotypes, (size_t) sizeof(char) );
		
			if(!visited_genotypes)fprintf(ferr, "count_genotypes: cannot c-allocate visited_genotypes, bye\n"), exit(3);

		
			if (opt_to== -1 && opt_from!= -1) // draw the landscape from opt_from
				{
				link_fitter_neighbor( f, fl, coord, width_geno, threshold, opt_from, geno_ref, R, opt_clear ,visited_genotypes,flat );
	//		exit (1);
				}
			else	
				if (opt_to!= -1 && opt_from == -1)   // draw the landscape to opt_to
					link_cripple_neighbor( f, fl, coord, width_geno, threshold, opt_to, geno_ref, R, opt_clear  ,visited_genotypes,flat);
		 		
				else{
				
				char *visited_genotypes = (char *)calloc( (size_t) fl->ngenotypes, (size_t) sizeof(char) );
				visited_genotypes[opt_to]=2;
//				link_from_to_neighbor( f, fl, coord, width_geno, threshold,opt_from,opt_to, geno_ref, R ,0 ,chemin, opt_clear,visited_links);
				link_from_to_neighbor( f, fl, coord, width_geno, threshold,opt_from,opt_to, geno_ref, R ,0 ,chemin,opt_clear,visited_genotypes,flat);
			
				free(visited_genotypes);
				}
			free(visited_genotypes);
		}
	/*
		Draw genotypes in the rotated plan
	*/
	

		for (g=fl->ngenotypes-1;g>=0;g--)
	{
//		genotype = int2genotype( *fl, g, genotype);

		if (fl->fitness[g] != DEFAULT_FITNESS){

			draw_geno(web,f, fl, g, coord, opt_log,R ,opt_from,opt_to,compact,flat);
		}
	}

	
	
	
	/*
		Close file
	*/
	if (with_legend==1)
	{
		int xgeno=70,
			xfrom=150,
			xto=200,
			
			xpeak=250,
			xsink=300,

			xhigh=350,
			xlow=400,

	
			xlinks=450,
	
	
		x=width;
		y=height+(2*margin_h)+45;
		
		SVG_draw_text(f,10,y,"Legend:",7);
		
		int sizeoftriange=10;
		y+=10;
		if (!compact)
			{
			for (i=0;i<3;i++)
//			fprintf(f,"<circle id=\"leg%d\" cx=\"%d\" cy=\"%d\"  r=\"%d\"  style=\" fill:#c6c7c8\" />\n",
//					i,100+(i*2*R)+1,y,2*R/3);
			fprintf(f,"<circle id=\"leg%d\" cx=\"%d\" cy=\"%d\"  r=\"%d\"  style=\" fill:#c6c7c8\" />\n",
					i,xgeno+(i*2*R)+1,y,2*R/3);

//			fprintf(f,"<circle id=\"leg%d\" cx=\"%d\" cy=\"%d\"  r=\"%d\"  style=\" fill:#1b171b\" />\n",
//					fl->nlocus,100+(i*2*R)+1,y,2*R/3);
			fprintf(f,"<circle id=\"leg%d\" cx=\"%d\" cy=\"%d\"  r=\"%d\"  style=\" fill:#1b171b\" />\n",
					4,xgeno+(i*2*R)+1,y,2*R/3);
			}		
		else
//			 fprintf(f,"<rect   x=\"%d\" height=\"15\" width=\"30\" y=\"%d\" rx=\"3\" ry=\"3\" style=\"fill:#cccccc;\" />\n", 
//		           100,y );
			 fprintf(f,"<rect   x=\"%d\" height=\"15\" width=\"30\" y=\"%d\" rx=\"3\" ry=\"3\" style=\"fill:#cccccc;\" />\n", 
		           xgeno,y );
		           
		SVG_draw_text(f,xgeno,y+40,"Genotype",5);
		
		SVG_draw_text(f,xsink,y+40,"Sink",5);
		SVG_draw_text(f,xpeak,y+40,"Peak",5);
		if (flat)
		{
		fprintf(f,"<rect id=\"leg12\"   x=\"%d\" y=\"%d\" height=\"15\" width=\"30\"  rx=\"3\" ry=\"3\" style=\"fill:#F8F8F8;stroke:#59ac2a;stroke-width:3;\" />\n", 
		          xpeak,y);
		fprintf(f,"<rect id=\"leg12\"   x=\"%d\" y=\"%d\" height=\"15\" width=\"30\"  rx=\"3\" ry=\"3\" style=\"fill:#181818;stroke:#e0141c;stroke-width:3;\" />\n", 
		           xsink,y);
	
		}
		else
		{
		 fprintf(f,"<rect id=\"leg12\"   x=\"%d\" y=\"%d\" height=\"15\" width=\"30\"  rx=\"3\" ry=\"3\" style=\"fill:#59ac2a;\" />\n", 
		           xpeak,y);
		 fprintf(f,"<rect id=\"leg13\"   x=\"%d\" y=\"%d\" height=\"15\" width=\"30\" rx=\"3\" ry=\"3\" style=\"fill:#e0141c;\" />\n", 
		           xsink,y);
		}
		if(flat)
		{
		SVG_draw_text(f,xhigh,y+40,"High",5);
		SVG_draw_text(f,xlow,y+40,"Low",5);

		fprintf(f,"<rect id=\"leg13\"   x=\"%d\" y=\"%d\" height=\"15\" width=\"30\" rx=\"3\" ry=\"3\" style=\"fill:#1B1740;\" />\n", 
		           xlow,y);
		fprintf(f,"<rect id=\"leg12\"   x=\"%d\" y=\"%d\" height=\"15\" width=\"30\"  rx=\"3\" ry=\"3\" style=\"fill:#F8F8F8;stroke:#181818;stroke-width:1;\" />\n", 
		           xhigh,y);           
		}
		if (opt_from!=-1)
		{
		SVG_draw_text(f,xfrom,y+40,"From",5);
		fprintf(f,"<path d=\"M %d %d l %d,-%d l 0,+%d l -%d,-%d \" stroke = \"red\" stroke-width = \"1\" fill = \"red\"/>",
				   xfrom,y,sizeoftriange,sizeoftriange,2*sizeoftriange,sizeoftriange,sizeoftriange);
		}	   
		if (opt_to!=-1)	
		{
			SVG_draw_text(f,xto,y+40,"To",5);
			fprintf(f,"<path d=\"M %d %d l %d,-%d l 0,+%d l -%d,-%d \" stroke = \"green\" stroke-width = \"1\" fill = \"green\"/>",
				   xto,y,sizeoftriange,sizeoftriange,2*sizeoftriange,sizeoftriange,sizeoftriange);
		}

	
		for (i=0;i<=2;i++)
		{

			SVG_draw_line(f,xlinks,y+(i*15),xlinks+30,y+(i*15),color[i],1);
			SVG_draw_text(f,xlinks+30+5,y+(i*15),legend[i],6);
		}
	}
/*	for (g=0;g<fl->nlocus;g++)
		{
		char leg[12];
		x= ((g+1) * width_geno * 3);
		SVG_draw_dashline(f,x,10,x,heighthg,"grey");
		}
*/
	if (web==0 && *(stats) != '\0')/*save the stats with the drawing*/
		{
		int xl;
		if (width+10 + 2*margin_w >450)
			xl=width+10 + 2*margin_w;
		else
			xl=570;	
		fprintf(f,"<text  x=\"%d\" y=\"20\" style=\"font-size:10pt;font-weight:bold;\">Statistiques</text>\n",xl );

		html_to_SVG(f,stats,xl,40, wstats);
		
//		html_to_SVG(f,stats,width+add_width+ 2*margin_w +wstats,40, wstats);
		}
	fprintf(f,"</g> </svg>\n");
	if (web)
	fprintf(f,"</div>\n");
	else
	fclose (f);

	free( genotype );
	free( geno_ref );
	free (chemin);
//	
}




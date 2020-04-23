#include <stdio.h>
#include <stdlib.h>
#define R 10 // rayon

long int fact(n)
{
if (n==1)
	return 1;
else
return (n*fact(n-1));
}

int max_trace(int nbloci)// only for 2 gives max of comb possible...
{
int i,j;
int n=0;
for (i=1;i<nbloci;i++)
	{
	j= fact(nbloci) /(float)(fact(i) * fact(nbloci-i));
	printf("%d\n",j);
 	if (j >n) n=j;
 	}
return n; 	
}


void svg_header(FILE *f,int w, int h,char *s)
{
fprintf(f,"<?xml version=\"1.0\"?>\n");
fprintf(f,"<svg width=\"%d\" height=\"%d\">\n",w,h);
fprintf(f,"<title>%s</title>\n",s);
}


void read (FILE *f)
{
	char *ligne,*ptr, *sep=" ,\t";
	int nbl=0;
	char **pheno;
	float *fitness;
	fscanf(f,"%[^\n]\n",ligne);
	ptr = strtok(ligne, sep);
	while ( ptr != NULL ) {
		 nbl++;
		 nba=atoi(ptr); //je dis quil ont tous le meme nbr dalles...
	     ptr = strtok(NULL, sep);
	}


	nbp=pow(nba,nbl);
	pheno=(char *)malloc(sizeof(char *) *nbp );
	for (i=0;i>nbp;i++)
		pheno[i]=malloc(sizeof(char)*nbl);
	for (i=0;i<nbp;i++)
	{
		for (j=0;j<nbl;j++)
			fscanf(f,"%d",&pheno[i][j]);
		fscanf(f,"%f",&fitness[i]);	
	}

}
/*draws a line of circles 1st one center is x,y*/
void draw_geno(FILE *f,int x, int y,int *pheno, float f, int nbloci)
{
	int i,j;
	char *colors[2]={"white","black"};
//char *colors[15]={"white","black","aqua","blue", "fuchsia", "gray", "green", "lime","maroon", "navy", "olive", "purple", "red", "silver", "teal", "yellow"};
	for (i=0;i<nbloci;i++)
	{
		fprintf(f,"<circle cx=\"%d\" cy=\"%d\" r=\"%d\" style=\"stroke: black; fill: %s\" />",x+(i*2*R)+1,y,R,colors[(int)pheno[i]]);		
	}
	fprintf(f,"<text x=\"%d\" y=\"%d\"style=\"text-anchor: middle\">%f</text>\n",x+((nbloci*R)/2.0),y+R+12,f);		
}

/*draws a set of genotype with same nubr of mut*/
}

void drawLand(struct landscape land,FILE *f)
{
int c,*g;
FILE *f;
int x,y;
int width,height,middle;
int width_geno=R*land.nlocus;
int height_geno=R*5;

int margin_left=10;
int margin_up=10;

width=width_geno*land.nlocus*2;
height=R*land.nlocus*land.nlocus*4;
middleY=height/2.0;
void svg_header(f,width, height,"what a beautifull view")

	for (i=0;i<land.ngenotypes;i++)
		{
		c=int2nbmut(land,i); 
		g=int2genotype( land, i, g);
		x=margin_left+(c*2*width_geno);
		y=margin_up+(c*2*height_geno);
		if (land.fitness[i]>=0)
		draw_geno(f,x, y,g, land.fitness[i],land.nlocus)

		}
fprintf(f,"</SVG>\n");

}


int main(int argc,char **argv)
{
	char *file, type;
	FILE *f;
	struct landscape land;
	
	if (argc!=3) printf("Usage %s [G|g] <nom fichier>\n",argv[0]);
	type =*argv[1]; 
	if (type !='g' && type !='G') printf("Usage %s [G|g] <nom fichier>\n",argv[0]);
	
	land= ReadFile (argv[2]);
	sprintf(file_out,"%s.svg",argv[2]);
	f=fopen(file_out,"w");
	if (type=='g')
		drawLand(land,f)
	else
		drawGraph(land,f);
		
	
	exit(0);
	
}


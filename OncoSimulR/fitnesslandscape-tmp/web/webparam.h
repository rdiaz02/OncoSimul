#ifndef _WEBPARAM_H_
#define _WEBPARAM_H_

typedef struct allparam
{
char *LAND_MODEL; // model name see HTML page to get all names or file name if adquate
int LAND_MODELNBR; //nbr of model
int LAND_NBLOCI; //nbr of gentotype
int LAND_CHAINS; //1 if only chains are drawn
int *ALLELS; //how many allels in each locus
char LAND_CHECK_SAMENBR;	// 0 if differents 1 if same

char 	LAND_LOG;// 0 is not log 1 if log
int 	LAND_REFERENCE; //nbr of left ref 
float 	LAND_THRESHOLD; //ratio fitnees
char	LAND_OPTCLEAN;// dsiplay neutrals
int 	LAND_DRAWFROM; //display lnks from
int 	LAND_DRAWTOEND; //display lnks to
float 	LAND_VALUE_MULTSAME; //case multiplicat
char 	LAND_MISSING; //1 if 0 
float 	LAND_SCALE_W; //width scale
float 	LAND_SCALE_H; //heigth scale
char 	*LAND_OUTFILE; 
char 	*LAND_DEFAULT;
char 	*LAND_PREVIOUS_FILE ;
int 	LAND_MULTSAME; //0 if differents value 1 othewise
int 	nbFields; // not exactly a parameter.. but how many fields decoded in the cgi input [usefull..]
float	LAND_RMFWEIGHT;
int 	LAND_WINDOW_WIDTHSIZE; //size of window
int 	LAND_WINDOW_HEIGHTSIZE; //size of window
int		LAND_ONLY_MUT;
int		LAND_COMPACT;
int		LAND_FLAT;
int 	*LAND_MASK;
} AllParams;


#endif

#ifndef _DRAWINGS_H_
#define _DRAWINGS_H_




int dist_to_genoref( struct landscape *h, int *g, int *ref );
int RetrieveRightNeighbors( struct landscape *fl, int *genotype, int *geno_ref, struct list *li );
void SVG_header(FILE *f,int w, int h,char *s);
void SVG_draw_line(FILE *f,int x_from,int y_from, int x_to,int y_to,char *color,int size);
void SVG_draw_dashline(FILE *f,int x1,int y1, int x2,int y2,char *color);
void SVG_draw_text(FILE *f,int x1,int y1,char *s,int R);
void SVG_draw_float(FILE *f,int x1,int y1,float v,int size, char *color);
void SVG_draw_int(FILE *f,int x1,int y1,int v,int size);
char *cold2warm( float frac, char *color);
//char *cold2warm( float frac);

void draw_geno(int web,FILE *f, struct landscape *fl, int g, int **coord, char opt_log,int R,int from,int to,int compact,int flat);
void draw_link(FILE *f, struct landscape *fl, int *geno1, int *geno2, int *geno_ref,  int **coord, int geno_width, char opt_neutral ,int R,float threshold,int flat);

void draw_links( FILE *f, struct landscape *fl, int **coord, int geno_width, float threshold, int opt_log, int opt_clear, int *geno_ref,int R, int OnlyMutation,int flat);
void link_fitter_neighbor( FILE *f, struct landscape *fl, int **coord, int geno_width,float threshold, int g, int *geno_ref  ,int R, int opt_clear, char *visited_genotypes,int flat);

void draw_FL( int web, struct landscape *fl,  char *file_out,  float threshold,  int opt_log,  int opt_clear, int opt_ref,int opt_to, int genoref, float rescale_height ,int R, int OnlyMutation,int with_legend,int opt_chains,int w,int h,int compact,int flat,char *);
void parseHTML(char *stats);
void html_to_SVG(FILE *f,char *stats,int x,int y, int w);
int trim(char *l);
#endif

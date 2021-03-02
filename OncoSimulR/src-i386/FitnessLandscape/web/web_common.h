#ifndef _WEBCOMMON_H_
#define _WEBCOMMON_H_



/*fns in web_common.c*/
void freeparam(AllParams *param);
void initparam(AllParams *param);
void save_SVG( struct landscape *fl,  char *file_out, AllParams *param,int R,char *stats);
void save_Landscape ( struct landscape *fl,  char *file_out);
void draw_menu (struct landscape land,  char *where, AllParams param, char *st,struct model_opt allopt);
void write_html_model( AllParams  param,struct model_opt allopt);
void draw_end();
void draw_header();
void setminMax(struct landscape *l);
void outputparam(AllParams param,struct model_opt myoptions);
void draw_lastline(AllParams param);
#endif

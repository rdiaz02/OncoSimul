
# ATTENTION pour pouvoir sauvegarder en pdf installer ca.. sudo apt-get install librsvg2-bin


BIN= fl_statistics fl_generate fl_draw  fl_genchains  fl_webdraw.cgi fl_webmodel.cgi fl_drawfile.cgi #test_gamma fl_gamma



HOME_BIN=~/bin/

# SGI:    RANLIB=touch  
# LINUX:  RANLIB=ranlib 
# MACOSX: RANLIB=ranlib 
RANLIB=ranlib

CC = gcc
ALLCFLAGS= -O -Wall 

SHELL=/bin/sh
LIBS = liblandscape.a
LIBS_DIR = -L.
LFLAGS = -llandscape

INCLUDE_DIR = -I.

AR = ar

SRC = gamma.c\
      landscape.c \
      random.c\
      sort.c\
      chain.c\
      input.c\
      common_drawings.c\
      genotypes.c\
      ordered_pairs.c\
      calculus.c\
      vector.c\
      generalized_chain.c\
      LinearAlgebra.c\
      summary_statistics.c\
      epistasis_type.c\
      decomposition.c\
      web/web_common.c\
      models.c

OBJECTS = $(SRC:.c=.o)


all: liblandscape.a $(BIN)

%.o: %.c
	$(CC) $(ALLCFLAGS) $(INCLUDE_DIR)  -c -o $@ $<;

#
# Check these two and move them into fl_real
# 

#test_gamma: test_gamma.c $(LIBS)
#	$(CC) $(ALLCFLAGS) -o $@ $(INCLUDE_DIR) $(LIBS_DIR) $< $(LFLAGS) -lm;

fl_genchains: main_generalizedchains_test.c $(LIBS)
	$(CC) $(ALLCFLAGS) -o $@ $(INCLUDE_DIR) $(LIBS_DIR) $< $(LFLAGS) -lm #-g

#
# Have a look and update if nescessary - very likelily deprecated
#

#fl_gamma: test_gamma.c $(LIBS)
#	$(CC) $(ALLCFLAGS) -o $@ $(INCLUDE_DIR) $(LIBS_DIR) $< $(LFLAGS) -lm


#
# Command line exe
# 

fl_statistics: main_fl_statistics.c $(LIBS)
	$(CC) $(ALLCFLAGS) -o $@ $(INCLUDE_DIR) $(LIBS_DIR) $< $(LFLAGS) -lm;

fl_generate: main_fl_generate.c  $(LIBS)
	$(CC) $(ALLCFLAGS) -o $@ $(INCLUDE_DIR) $(LIBS_DIR) $< $(LFLAGS) -lm;

fl_draw: main_draw.c $(LIBS)
	$(CC) $(ALLCFLAGS) -o $@ $(INCLUDE_DIR) $(LIBS_DIR) $< $(LFLAGS) -lm;


#
# Web exe
# 

fl_webdraw.cgi: main_webdraw.c  $(LIBS)
	$(CC) $(ALLCFLAGS) -o $@ $(INCLUDE_DIR) $(LIBS_DIR) $< $(LFLAGS) -lm;

fl_webmodel.cgi: main_webmodels.c   $(LIBS)
	$(CC) $(ALLCFLAGS) -o $@ $(INCLUDE_DIR) $(LIBS_DIR) $< $(LFLAGS) -lm;

fl_drawfile.cgi:	web/draw_custom.c  $(LIBS)
	$(CC) $(ALLCFLAGS) -o $@ $(INCLUDE_DIR) $(LIBS_DIR) $< $(LFLAGS) -lm;

#
# Librarys
#

liblandscape.a: $(OBJECTS)
	$(AR) -r $@ $+ ;
	$(RANLIB) $@;

install: all
	cp $(BIN) $(HOME_BIN)


installmac:
	@ if [ -f /usr/local/bin/rsvg-convert ]; then echo ""; else echo "**************";echo "**********!!!!! install rsvg-convert in /usr/local/bin/**********!!!!! ";echo "**************";fi
	zip fl_formac.zip fl_draw fl_statistics fl_generate fl_genchains
	cp fl_formac.zip /Applications/MAMP/htdocs/land/
	cp fl_webdraw.cgi /Applications/MAMP/cgi-bin/
	cp fl_webmodel.cgi /Applications/MAMP/cgi-bin/
	cp fl_drawfile.cgi /Applications/MAMP/cgi-bin/
	cp web/land_js.js /Applications/MAMP/htdocs/land/js
	cp web/landmodel_js.js /Applications/MAMP/htdocs/land/js
	cp web/land.css /Applications/MAMP/htdocs/land/
	cp web/magellan.css /Applications/MAMP/htdocs/land/
	cp web/Magellan.main.mac.html /Applications/MAMP/htdocs//land/
	cp web/Magellan.help.mac.html /Applications/MAMP/htdocs//land/
	cp web/Magellan.custom.mac.html /Applications/MAMP/htdocs/land/
	cp web/Magellan.model.mac.html /Applications/MAMP/htdocs/land/
	cp web/Magellan.known.mac.html /Applications/MAMP/htdocs/land/

	
installserver:
	@ if [ -f /usr/bin/rsvg-convert ]; then echo ""; else echo "**************";echo "**********!!!!! install rsvg-convert in /usr/bin/**********!!!!! ";echo "**************";fi
	zip fl_forunix.zip fl_draw fl_statistics fl_generate fl_genchains
	cp fl_forunix.zip /var/www/html/wwwabi/public/Magellan
	cp fl_webdraw.cgi /var/www/cgi-bin/
	cp fl_webmodel.cgi /var/www/cgi-bin/
	cp fl_drawfile.cgi /var/www/cgi-bin/
	cp web/land_js.js /var/www/html/wwwabi/public/Magellan/js
	cp web/landmodel_js.js /var/www/html/wwwabi/public/Magellan/js
	cp web/land.css /var/www/html/wwwabi/public/Magellan/
	cp web/magellan.css /var/www/html/wwwabi/public/Magellan/
	cp web/Magellan.main.html /var/www/html/wwwabi/public/Magellan
	cp web/Magellan.help.html /var/www/html/wwwabi/public/Magellan
	cp web/Magellan.custom.html /var/www/html/wwwabi/public/Magellan
	cp web/Magellan.model.html /var/www/html/wwwabi/public/Magellan
	cp web/Magellan.known.html /var/www/html/wwwabi/public/Magellan
	
archive: clean
	\rm -rf FitnessLandscape; mkdir FitnessLandscape;
	cp -r *.h *.c web Makefile FitnessLandscape ;
	tmp=`date | awk '{print $$2$$3$$4}'`; tar czvf ../FitnessLandscape.$$tmp.tgz FitnessLandscape;
	\rm -rf FitnessLandscape

smallarchive: clean
	\rm -rf FitnessLandscape; mkdir FitnessLandscape;
	cp *.h *.c  Makefile FitnessLandscape ;
	cp smallMake FitnessLandscape/Makefile
	tmp=`date | awk '{print $$2$$3$$4}'`; tar czvf ../MagellanCL.$$tmp.tgz FitnessLandscape;
	\rm -rf FitnessLandscape

clean:
	\rm -f *.o liblandscape.a $(BIN)
	\rm -f *.cgi
	\rm -rf *.dSYM
	\rm -f ./web/*.o
	\rm -f *~
	\rm -f GTAGS
	\rm -f GRTAGS
	\rm -f GPATH


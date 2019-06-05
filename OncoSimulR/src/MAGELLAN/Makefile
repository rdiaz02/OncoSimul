
# ATTENTION pour pouvoir sauvegarder en pdf installer ca.. sudo apt-get install librsvg2-bin


## BIN= fl_statistics fl_generate fl_draw  fl_genchains  fl_webdraw.cgi fl_webmodel.cgi fl_drawfile.cgi #test_gamma fl_gamma
## RDU
BIN= fl_statistics fl_generate fl_draw  fl_genchains ## fl_webdraw.cgi fl_webmodel.cgi fl_drawfile.cgi #test_gamma fl_gamma


# HOME_BIN=~/bin/

# SGI:    RANLIB=touch  
# LINUX:  RANLIB=ranlib 
# MACOSX: RANLIB=ranlib 
RANLIB=ranlib

CC = gcc
ALLCFLAGS= -O -Wall 

## RDU
## SHELL=/bin/sh
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
      models.c

# OBJECTS = gamma.o\
#       landscape.o \
#       random.o\
#       sort.o\
#       chain.o\
#       input.o\
#       common_drawings.o\
#       genotypes.o\
#       ordered_pairs.o\
#       calculus.o\
#       vector.o\
#       generalized_chain.o\
#       LinearAlgebra.o\
#       summary_statistics.o\
#       epistasis_type.o\
#       decomposition.o\
#       models.o
OBJECTS = $(SRC:.c=.o)

# ## RDU
# magellanall: liblandscape.a fl_statistics fl_generate fl_genchains
# ## all: liblandscape.a $(BIN)

%.o: %.c
	$(CC) $(ALLCFLAGS) $(INCLUDE_DIR)  -c -o $@ $<;

#
# Librarys
#

liblandscape.a: $(OBJECTS)
	$(AR) -r $@ $+ ;
	$(RANLIB) $@;


fl_genchains: main_generalizedchains_test.c $(LIBS)
	$(CC) $(ALLCFLAGS) -o $@ $(INCLUDE_DIR) $(LIBS_DIR) $< $(LFLAGS) -lm #-g

fl_statistics: main_fl_statistics.c $(LIBS)
	$(CC) $(ALLCFLAGS) -o $@ $(INCLUDE_DIR) $(LIBS_DIR) $< $(LFLAGS) -lm;

fl_generate: main_fl_generate.c  $(LIBS)
	$(CC) $(ALLCFLAGS) -o $@ $(INCLUDE_DIR) $(LIBS_DIR) $< $(LFLAGS) -lm;






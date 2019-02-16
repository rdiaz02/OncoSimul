## Look at p. 29 and ff. of http://dirk.eddelbuettel.com/papers/useR2010hpcTutorial.pdf

## I also remember kcachegrind was nice to see output.

## Old example, with some confussed notes. But maybe a way to get going?

## using 2nd approach in Paul Johnson's help
## http://pj.freefaculty.org/blog/?p=140

## in Lacerta or Gallotia
R-2.15 CMD INSTALL MatherR -l /usr/local/sources/R-2.15.3-patched-03-03/library

/usr/bin/R CMD INSTALL MatherR -l /home/ramon/R/library

## export LD_PRELOAD=/usr/lib/libprofiler.so

CPUPROFILE="myprof.log" /usr/bin/R --no-save --no-restore < debug-call.R > debug-call.Rout

## R --no-save --no-restore < debug-call.R
google-pprof --text /usr/lib/R/bin/exec/R myprof.log  

google-pprof --text --focus=Algorithm5F /usr/lib/R/bin/exec/R myprof.log | more
google-pprof --text --focus=Algorithm5G /usr/lib/R/bin/exec/R myprof.log | more

google-pprof --text --focus=Algorithm5G --lines /usr/lib/R/bin/exec/R myprof.log | more



## get rid of the unkown stuff
google-pprof --text --focus=Algorithm5G /usr/lib/R/bin/exec/R myprof.log | grep -v 0000

## sort by samples in function
google-pprof --text --focus=Algorithm5K /usr/local/sources/R-2.15.3-patched-03-03/bin/exec/R myprof.log | grep -v 0000 | sort -k 5 -n


## google-pprof --gv --focus=Algorithm5D /usr/lib/R/bin/exec/R myprof.log 

google-pprof --callgrind --focus=Algorithm5G /usr/lib/R/bin/exec/R myprof.log > mp.callgrind; kcachegrind mp.callgrind


google-pprof --gv --focus=Algorithm5G /usr/lib/R/bin/exec/R myprof.log


## cache misses and memory stuff using valgrind
R -d "valgrind --tool=cachegrind" --no-save --no-restore < a1.R

R -d "valgrind --tool=callgrind --cache-sim=yes --branch-sim=yes" --no-save --no-restore < a1.R

## the following can fail 
R -d "valgrind --tool=massif" --no-save --no-restore < a1.R
## then ms_print and name of massif.out



#############  Might need to compile with specific Makevars, either nder
#############  .R or under the /src of the package

## but what follows is for profile usage/generation

## generate
#MYAMD="-march=native -O3 -g -fpic -pipe"
OPTIMFLAGS=-march=native -O3 -ffunction-sections -fprofile-generate -fprofile-arcs -g -fpic -pipe
CFLAGS=-march=native -O3 -ffunction-sections  -fprofile-generate -fprofile-arcs  -g -fpic -pipe
CXXFLAGS=-march=native -O3 -ffunction-sections  -fprofile-generate -fprofile-arcs -g -fpic -pipe
FFLAGS=-march=native -O3 -ffunction-sections  -fprofile-generate -fprofile-arcs -g -fpic -pipe
FCFLAGS=-march=native -O3 -ffunction-sections  -fprofile-generate -fprofile-arcs -g -fpic -pipe


### use
#MYAMD="-march=native -O3 -g -fpic -pipe"
OPTIMFLAGS=-march=native -O3 -ffunction-sections -fprofile-use -fprofile-arcs -g -fpic -pipe
CFLAGS=-march=native -O3 -ffunction-sections  -fprofile-use -fprofile-arcs  -g -fpic -pipe
CXXFLAGS=-march=native -O3 -ffunction-sections  -fprofile-use -fprofile-arcs -g -fpic -pipe
FFLAGS=-march=native -O3 -ffunction-sections  -fprofile-use -fprofile-arcs -g -fpic -pipe
FCFLAGS=-march=native -O3 -ffunction-sections  -fprofile-use -fprofile-arcs -g -fpic -pipe


#### In the package? But I think most unneded. Except add, in Makevars, -lgmpxx -lgmp -lgcov?
# combine with standard arguments for R
PKG_CPPFLAGS = `${R_HOME}/bin/Rscript -e "RcppGSL:::CFlags()"` -I../inst/include -Wall

## next is for profling with gcc
PKG_LIBS = `$(R_HOME)/bin/Rscript -e "Rcpp:::LdFlags()"` ` ${R_HOME}/bin/Rscript -e "RcppGSL:::LdFlags()"` $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) -lgmpxx -lgmp -lgcov

#!/bin/bash

V_R=$1

V_ADA=$(cat Full/OncoSimulR/DESCRIPTION | grep Version | cut -d' ' -f2)

rm OncoSimulR_$V_ADA.tar.gz

rm ./Full/OncoSimulR/src/*.so
rm ./Full/OncoSimulR/src/*.o
rm ./Full/OncoSimulR/src/OncoSimulR.so
rm ./Full/OncoSimulR/src/OncoSimul.o
rm ./Full/OncoSimulR/src/OncoSimulR_init.o
rm ./Full/OncoSimulR/src/symbols.rds
rm ./Full/OncoSimulR.Rcheck/* -r -f
rm ./OncoSimulR.Rcheck/* -r -f
# rm ./Full/OncoSimulR/inst/doc/auto/*
# rmdir ./Full/OncoSimulR/inst/doc/auto
rm ./Full/OncoSimulR/vignettes/auto/*
rmdir ./Full/OncoSimulR/vignettes/auto
rm ./Full/OncoSimulR/vignettes/*.bbl
rm ./Full/OncoSimulR/vignettes/*.aux
rm ./Full/OncoSimulR/vignettes/*.toc
rm ./Full/OncoSimulR/vignettes/*.tex
rm ./Full/OncoSimulR/vignettes/*.pdf
rm ./Full/OncoSimulR/vignettes/*.log
rm ./Full/OncoSimulR/vignettes/*.out
rm ./Full/OncoSimulR/vignettes/*.blg



time $V_R CMD build --keep-empty-dirs --resave-data Full/OncoSimulR

time $V_R CMD check --no-vignettes --timings OncoSimulR_$V_ADA.tar.gz

time $V_R CMD check --as-cran --no-vignettes --timings OncoSimulR_$V_ADA.tar.gz


## --force-multiarch in windoze?

#!/bin/bash

V_R=$1

V_ADA=$(cat ./OncoSimulR/DESCRIPTION | grep Version | cut -d' ' -f2)

rm OncoSimulR_$V_ADA.tar.gz

rm ./OncoSimulR/src/*.so
rm ./OncoSimulR/src/*.o
rm ./OncoSimulR/src/OncoSimulR.so
rm ./OncoSimulR/src/OncoSimul.o
rm ./OncoSimulR/src/OncoSimulR_init.o
rm ./OncoSimulR/src/symbols.rds
rm /OncoSimulR.Rcheck/* -r -f
rm ./OncoSimulR.Rcheck/* -r -f
# rm ./OncoSimulR/inst/doc/auto/*
# rmdir ./OncoSimulR/inst/doc/auto
rm ./OncoSimulR/vignettes/auto/*
rmdir ./OncoSimulR/vignettes/auto
rm ./OncoSimulR/vignettes/figure/*
rmdir ./OncoSimulR/vignettes/figure
rm ./OncoSimulR/vignettes/*.bbl
rm ./OncoSimulR/vignettes/*.aux
rm ./OncoSimulR/vignettes/*.toc
rm ./OncoSimulR/vignettes/*.tex
rm ./OncoSimulR/vignettes/*.pdf
rm ./OncoSimulR/vignettes/*.log
rm ./OncoSimulR/vignettes/*.out
rm ./OncoSimulR/vignettes/*.blg
rm ./OncoSimulR/vignettes/*.synctex.*


echo "\n **************** \n"
echo "************ R CMD build\n"
time $V_R CMD build --keep-empty-dirs OncoSimulR
## time $V_R CMD build --keep-empty-dirs --resave-data OncoSimulR
echo "\n **************** \n"
echo "R CMD check\n"
time $V_R CMD check --no-vignettes --timings OncoSimulR_$V_ADA.tar.gz
## time $V_R CMD check --as-cran --no-vignettes --timings OncoSimulR_$V_ADA.tar.gz
echo "\n **************** \n"
echo "long manual tests\n"
$V_R CMD INSTALL OncoSimulR_$V_ADA.tar.gz
time $V_R -e 'library(OncoSimulR); library(testthat); test_dir("./OncoSimulR/tests/manual/")'

## --force-multiarch in windoze?

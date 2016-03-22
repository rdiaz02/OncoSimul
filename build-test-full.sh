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


echo " ************************************** "
echo " **********   R CMD build   *********** "
echo ""
time $V_R CMD build --keep-empty-dirs OncoSimulR
echo " "
echo " =======      done R CMD build   ======= "
echo " "
## time $V_R CMD build --keep-empty-dirs --resave-data OncoSimulR
echo " ************************************** "
echo " *********     R CMD check   ********** "
echo " "
time $V_R CMD check --no-vignettes --timings OncoSimulR_$V_ADA.tar.gz
echo " "
echo " =======      done R CMD check   =======  "
echo " "
## time $V_R CMD check --as-cran --no-vignettes --timings OncoSimulR_$V_ADA.tar.gz
echo " ******************************************* "
echo " *******      long manual tests ************ "
$V_R CMD INSTALL OncoSimulR_$V_ADA.tar.gz
time $V_R -e 'library(OncoSimulR); library(testthat); test_dir("./OncoSimulR/tests/manual/")'
echo " "
echo " =======     done long manual tests   =======     "
echo " "
## --force-multiarch in windoze?

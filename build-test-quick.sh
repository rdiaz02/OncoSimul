#!/bin/bash

## Do not build vignettes. Allows quick tests 
V_R=$1

V_ADA=$(cat ./OncoSimulR/DESCRIPTION | grep Version | cut -d' ' -f2)

rm OncoSimulR_$V_ADA.tar.gz

rm ./OncoSimulR/Rplots.pdf
rm ./OncoSimulR/.Rhistory
rm ./OncoSimulR/*~
rm ./OncoSimulR/inst/*~
rm ./OncoSimulR/man/*~
rm ./OncoSimulR/man/Rplots.pdf
rm ./OncoSimulR/man/.Rhistory
rm ./OncoSimulR/tests/*~
rm ./OncoSimulR/tests/Rplots.pdf
rm ./OncoSimulR/tests/.Rhistory
rm ./OncoSimulR/data/*~
rm ./OncoSimulR/data/Rplots.pdf
rm ./OncoSimulR/data/.Rhistory
rm ./OncoSimulR/tests/testthat/*~
rm ./OncoSimulR/tests/testthat/Rplots.pdf
rm ./OncoSimulR/tests/testthat/.Rhistory
rm ./OncoSimulR/tests/manual/*~
rm ./OncoSimulR/tests/manual/Rplots.pdf
rm ./OncoSimulR/tests/manual/.Rhistory
rm ./OncoSimulR/R/*~
rm ./OncoSimulR/R/.Rhistory
rm ./OncoSimulR/R/Rplots.pdf
rm ./OncoSimulR/vignettes/*~
rm ./OncoSimulR/vignettes/Rplots.pdf
rm ./OncoSimulR/vignettes/.Rhistory
rm ./OncoSimulR/src/*.so
rm ./OncoSimulR/src/*~
rm ./OncoSimulR/src/*.o
rm ./OncoSimulR/src/*.gcda
rm ./OncoSimulR/src/*.gcno
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
rm ./OncoSimulR/vignettes/OncoSimulR.tex
rm ./OncoSimulR/vignettes/*.pdf
rm ./OncoSimulR/vignettes/*.log
rm ./OncoSimulR/vignettes/*.out
rm ./OncoSimulR/vignettes/*.blg
rm ./OncoSimulR/vignettes/*.synctex.*
rm ./OncoSimulR/src/liblandscape.a
rm ./OncoSimulR/src/fl_statistics 
rm ./OncoSimulR/src/fl_generate
rm ./OncoSimulR/src/fl_genchains
rm ./OncoSimulR/src/FitnessLandscape/*.o
rm ./OncoSimulR/src/FitnessLandscape/*~
make -C ./OncoSimulR/src/FitnessLandscape clean

## We always do this, though it should not be necessary
sed -i 's/^};$/}/' ./OncoSimulR/src/FitnessLandscape/generalized_chain.c 

export R_CHECK_ENVIRON="~/.R/check.Renviron"
echo " ************************************** "
echo " **********   R CMD build --vanilla  *********** "
echo ""
time $V_R --vanilla CMD build --no-build-vignettes --no-manual --no-resave-data --keep-empty-dirs OncoSimulR
echo " "
echo " =======      done R CMD build --vanilla  ======= "
echo " "
echo " ************************************** "
echo " *********     R CMD check --vanilla  ********** "
echo " "
time $V_R --vanilla CMD check --no-vignettes --no-build-vignettes --ignore-vignettes --no-manual --timings OncoSimulR_$V_ADA.tar.gz
echo " "
echo " =======      done R CMD check --vanilla  =======  "
echo " "
echo " =======  installing with tests "
$V_R --vanilla  CMD INSTALL --install-tests OncoSimulR_$V_ADA.tar.gz

    

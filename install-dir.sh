#!/bin/bash

V_R=$1
# export R_MAKEVARS_USER=/home/ramon/.R/Makevars
# export CXX="clang++"
# export CC="clang"
# export CXX1X="clang++"
# export CFLAGS=" -Wall -g -O2"
# export CPICFLAGS=" -fPIC"
# export CPP="clang -E"
# export CXXFLAGS=" -Wall -g -O2"
# export CXXPICFLAGS=" -fPIC"
# export CXXCPP="clang++ -E"
# export CXX1XFLAGS=" -Wall -g -O2"
# export CXX1XPICFLAGS=" -fPIC"
# export CXX1XSTD=" -std=gnu++11"

rm ./OncoSimulR/src/OncoSimulR.so
rm ./OncoSimulR/src/OncoSimul.o
rm ./OncoSimulR/src/OncoSimulR_init.o
rm ./OncoSimulR/src/symbols.rds
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
rm ./OncoSimulR/src/OncoSimulR.so
rm ./OncoSimulR/src/OncoSimul.o
rm ./OncoSimulR/src/OncoSimulR_init.o
rm ./OncoSimulR/src/symbols.rds
rm ./OncoSimulR/src/liblandscape.a
rm ./OncoSimulR/src/fl_statistics 
rm ./OncoSimulR/src/fl_generate
rm ./OncoSimulR/src/fl_genchains
rm ./OncoSimulR/src/FitnessLandscape/liblandscape.a
rm ./OncoSimulR/src/FitnessLandscape/*.o
rm ./OncoSimulR/src/FitnessLandscape/*~
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
# rm ./OncoSimulR/vignettes/*.tex  ## NO!
rm ./OncoSimulR/vignettes/*.pdf
rm ./OncoSimulR/vignettes/*.log
rm ./OncoSimulR/vignettes/*.out
rm ./OncoSimulR/vignettes/*.blg
rm ./OncoSimulR/vignettes/*.synctex.*

make -C ./OncoSimulR/src/FitnessLandscape clean

$V_R --vanilla CMD INSTALL --install-tests ./OncoSimulR


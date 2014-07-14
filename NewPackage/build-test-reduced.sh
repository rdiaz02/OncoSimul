#!/bin/bash

V_R=$1

V_ADA=$(cat Reduced/OncoSimulR/DESCRIPTION | grep Version | cut -d' ' -f2)

rm OncoSimulR_$V_ADA.tar.gz

rm ./Reduced/OncoSimulR/src/*.so
rm ./Reduced/OncoSimulR/src/*.o
rm ./Reduced/OncoSimulR/src/OncoSimulR.so
rm ./Reduced/OncoSimulR/src/OncoSimul.o
rm ./Reduced/OncoSimulR/src/OncoSimulR_init.o
rm ./Reduced/OncoSimulR/src/symbols.rds
rm ./Reduced/OncoSimulR.Rcheck/* -r -f
rm ./OncoSimulR.Rcheck/* -r -f
# rm ./Reduced/OncoSimulR/inst/doc/auto/*
# rmdir ./Reduced/OncoSimulR/inst/doc/auto
rm ./Reduced/OncoSimulR/vignettes/auto/*
rmdir ./Reduced/OncoSimulR/vignettes/auto
rm ./Reduced/OncoSimulR/vignettes/*.bbl
rm ./Reduced/OncoSimulR/vignettes/*.aux
rm ./Reduced/OncoSimulR/vignettes/*.toc
rm ./Reduced/OncoSimulR/vignettes/*.tex
rm ./Reduced/OncoSimulR/vignettes/*.pdf
rm ./Reduced/OncoSimulR/vignettes/*.log
rm ./Reduced/OncoSimulR/vignettes/*.out
rm ./Reduced/OncoSimulR/vignettes/*.blg



time $V_R CMD build --keep-empty-dirs --no-resave-data Reduced/OncoSimulR

time $V_R CMD check --no-vignettes --timings OncoSimulR_$V_ADA.tar.gz

time $V_R CMD check --as-cran --no-vignettes --timings OncoSimulR_$V_ADA.tar.gz


## --force-multiarch in windoze?

## can use BiocCheck: install, load, and run with dir
## https://github.com/Bioconductor/BiocCheck
## BiocCheck("~/Proyectos/OncoSimulR/NewPackage/Reduced/OncoSimulR/")
## I hack it to provide more detailed output SOurces/BiocCheckR, as library BiocCheckRamon

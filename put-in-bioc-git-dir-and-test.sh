#!/bin/bash

## This places the code to upload to BioC git repo in the right place.

## If passed an additional argument, with the path of an R version, it
## builds and tests.

## Yes, works under eshell.


## I need to clean up first, because otherwise when I remove files in the
## master git repo, they are not necessarily removed in the svn repo

rm ../BioConductor-git/OncoSimulR/vignettes/*
rm ../BioConductor-git/OncoSimulR/tests/*
rm ../BioConductor-git/OncoSimulR/tests/manual/*
rm ../BioConductor-git/OncoSimulR/tests/testthat/*
rm ../BioConductor-git/OncoSimulR/src/*
rm ../BioConductor-git/OncoSimulR/R/*
rm ../BioConductor-git/OncoSimulR/man/*
rm ../BioConductor-git/OncoSimulR/data/*
rm ../BioConductor-git/OncoSimulR/inst/*
rm ../BioConductor-git/OncoSimulR/inst/miscell/*

cp OncoSimulR/vignettes/relfunct.tex ../BioConductor-git/OncoSimulR/vignettes/.
cp OncoSimulR/vignettes/preamble.tex ../BioConductor-git/OncoSimulR/vignettes/.
cp OncoSimulR/vignettes/relfunct.png ../BioConductor-git/OncoSimulR/vignettes/.
cp OncoSimulR/vignettes/custom4.css ../BioConductor-git/OncoSimulR/vignettes/.
cp OncoSimulR/vignettes/OncoSimulR.Rmd ../BioConductor-git/OncoSimulR/vignettes/.
cp OncoSimulR/vignettes/OncoSimulR.bib ../BioConductor-git/OncoSimulR/vignettes/.
cp OncoSimulR/vignettes/gitsetinfo.sty ../BioConductor-git/OncoSimulR/vignettes/.
cp OncoSimulR/vignettes/gitinfo.sty ../BioConductor-git/OncoSimulR/vignettes/.
cp OncoSimulR/vignettes/gitHeadInfo.gin ../BioConductor-git/OncoSimulR/vignettes/.

cp OncoSimulR/src/*.c ../BioConductor-git/OncoSimulR/src/.
cp OncoSimulR/src/*.cpp ../BioConductor-git/OncoSimulR/src/.
cp OncoSimulR/src/*.h ../BioConductor-git/OncoSimulR/src/.
cp OncoSimulR/src/Makevars* ../BioConductor-git/OncoSimulR/src/.

cp OncoSimulR/R/*.R ../BioConductor-git/OncoSimulR/R/.

cp OncoSimulR/tests/testthat.R ../BioConductor-git/OncoSimulR/tests/.
cp OncoSimulR/tests/testthat/*.R ../BioConductor-git/OncoSimulR/tests/testthat/.
cp OncoSimulR/tests/manual/*.R ../BioConductor-git/OncoSimulR/tests/manual/.
cp OncoSimulR/tests/manual/*.txt ../BioConductor-git/OncoSimulR/tests/manual/.
cp OncoSimulR/tests/*.txt ../BioConductor-git/OncoSimulR/tests/.

cp OncoSimulR/man/*.Rd ../BioConductor-git/OncoSimulR/man/.

cp OncoSimulR/inst/NEWS ../BioConductor-git/OncoSimulR/inst/.
cp OncoSimulR/inst/CITATION ../BioConductor-git/OncoSimulR/inst/.
cp OncoSimulR/inst/miscell/example-binom-problems.cpp ../BioConductor-git/OncoSimulR/inst/miscell/.
cp OncoSimulR/inst/miscell/*.R ../BioConductor-git/OncoSimulR/inst/miscell/.
cp OncoSimulR/inst/testdata_fee.RData ../BioConductor-git/OncoSimulR/inst/.

cp OncoSimulR/data/*.RData ../BioConductor-git/OncoSimulR/data/.

cp OncoSimulR/NAMESPACE ../BioConductor-git/OncoSimulR/.

cp OncoSimulR/DESCRIPTION ../BioConductor-git/OncoSimulR/.


## should we run the tests?
export R_CHECK_ENVIRON="~/.R/check.Renviron"
if [[ $# == 1 ]]; then
    V_R=$1
    cd ~/Proyectos/BioConductor-git
    V_P=$(cat ./OncoSimulR/DESCRIPTION | grep Version | cut -d' ' -f2)
    rm OncoSimulR_$V_P.tar.gz
    ## As shown in build report from BioC
    # echo " ***************************************** "
    # echo " *********  R CMD build   ************** "
    # echo " "
    # time $V_R CMD build --keep-empty-dirs --no-resave-data OncoSimulR
    # echo " "
    # echo " ===========  done R CMD build  ========== "
    # echo " "
    # ## As shown in check report from BioC
    # echo " *************************************************** "
    # echo " **** R CMD check , as in check report **** "
    # echo ""
    # time $V_R CMD check --no-vignettes --timings OncoSimulR_$V_P.tar.gz
    # echo " "
    # echo " =========   done R CMD check as in check report  =======  "
    # echo " "
    ## time as explained in https://www.bioconductor.org/developers/package-guidelines/#correctness
    # echo " ************************************ "
    # echo " *****   R CMD check: time OK?  ***** "
    # echo ""
    # time $V_R CMD check --no-build-vignettes OncoSimulR_$V_P.tar.gz
    # echo " "
    # echo " ===========  done R CMD check time OK?   ========== "
    ## As shown in build report from BioC
    echo " ***************************************** "
    echo " *********  R CMD build vanilla  ************** "
    echo " "
    time $V_R --vanilla CMD build --keep-empty-dirs --no-resave-data OncoSimulR
    echo " "
    echo " ===========  done R CMD build   vanilla ========== "
    echo " "
    ## As shown in check report from BioC
    echo " *************************************************** "
    echo " **** R CMD check vanilla , as in check report **** "
    echo ""
    time $V_R --vanilla CMD check --no-vignettes --timings OncoSimulR_$V_P.tar.gz
    echo " "
    echo " =========   done R CMD check vanilla as in check report  =======  "
    echo " "
    ## time as explained in https://www.bioconductor.org/developers/package-guidelines/#correctness
    echo " ************************************ "
    # echo " *****   R CMD check: time OK?  ***** "
    # echo ""
    # time $V_R --vanilla CMD check --no-build-vignettes OncoSimulR_$V_P.tar.gz
    # echo " "
    # echo " ===========  done R CMD check vanilla time OK?   ========== "

    
fi


## Check what/if things need adding
cd ~/Proyectos/BioConductor-git/OncoSimulR

git status | less


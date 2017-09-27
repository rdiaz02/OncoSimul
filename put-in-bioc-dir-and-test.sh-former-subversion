#!/bin/bash

## This places the code to upload to BiC in the right place.
## The Subversion-in-BioC directory is a special directory (under svn vc)
## for the code in BioC. Yes, will eventually use
## the bridge or similar.

## If passed an additional argument, with the path of an R version, it
## builds and tests.

## Yes, works under eshell.

## I need to clean up first, because otherwise when I remove files in the
## master git repo, they are not necessarily removed in the svn repo

rm ./Subversion-in-BioC/OncoSimulR/vignettes/*
rm ./Subversion-in-BioC/OncoSimulR/tests/*
rm ./Subversion-in-BioC/OncoSimulR/tests/manual/*
rm ./Subversion-in-BioC/OncoSimulR/tests/testthat/*
rm ./Subversion-in-BioC/OncoSimulR/src/*
rm ./Subversion-in-BioC/OncoSimulR/R/*
rm ./Subversion-in-BioC/OncoSimulR/man/*
rm ./Subversion-in-BioC/OncoSimulR/inst/*
rm ./Subversion-in-BioC/OncoSimulR/data/*
rm ./Subversion-in-BioC/OncoSimulR/inst/miscell/*

cp OncoSimulR/vignettes/relfunct.tex ./Subversion-in-BioC/OncoSimulR/vignettes/.
cp OncoSimulR/vignettes/preamble.tex ./Subversion-in-BioC/OncoSimulR/vignettes/.
cp OncoSimulR/vignettes/relfunct.png ./Subversion-in-BioC/OncoSimulR/vignettes/.
cp OncoSimulR/vignettes/custom4.css ./Subversion-in-BioC/OncoSimulR/vignettes/.
cp OncoSimulR/vignettes/OncoSimulR.Rmd ./Subversion-in-BioC/OncoSimulR/vignettes/.
cp OncoSimulR/vignettes/OncoSimulR.bib ./Subversion-in-BioC/OncoSimulR/vignettes/.
cp OncoSimulR/vignettes/gitsetinfo.sty ./Subversion-in-BioC/OncoSimulR/vignettes/.
cp OncoSimulR/vignettes/gitinfo.sty ./Subversion-in-BioC/OncoSimulR/vignettes/.
cp OncoSimulR/vignettes/gitHeadInfo.gin ./Subversion-in-BioC/OncoSimulR/vignettes/.

cp OncoSimulR/src/*.c ./Subversion-in-BioC/OncoSimulR/src/.
cp OncoSimulR/src/*.cpp ./Subversion-in-BioC/OncoSimulR/src/.
cp OncoSimulR/src/*.h ./Subversion-in-BioC/OncoSimulR/src/.
cp OncoSimulR/src/Makevars* ./Subversion-in-BioC/OncoSimulR/src/.

cp OncoSimulR/R/*.R ./Subversion-in-BioC/OncoSimulR/R/.

cp OncoSimulR/tests/testthat.R ./Subversion-in-BioC/OncoSimulR/tests/.
cp OncoSimulR/tests/testthat/*.R ./Subversion-in-BioC/OncoSimulR/tests/testthat/.
cp OncoSimulR/tests/manual/*.R ./Subversion-in-BioC/OncoSimulR/tests/manual/.
cp OncoSimulR/tests/manual/*.txt ./Subversion-in-BioC/OncoSimulR/tests/manual/.
cp OncoSimulR/tests/*.txt ./Subversion-in-BioC/OncoSimulR/tests/.

cp OncoSimulR/man/*.Rd ./Subversion-in-BioC/OncoSimulR/man/.

cp OncoSimulR/inst/NEWS ./Subversion-in-BioC/OncoSimulR/inst/.
cp OncoSimulR/inst/CITATION ./Subversion-in-BioC/OncoSimulR/inst/.
cp OncoSimulR/inst/miscell/example-binom-problems.cpp ./Subversion-in-BioC/OncoSimulR/inst/miscell/.
cp OncoSimulR/inst/miscell/*.R ./Subversion-in-BioC/OncoSimulR/inst/miscell/.

cp OncoSimulR/data/*.RData ./Subversion-in-BioC/OncoSimulR/data/.

cp OncoSimulR/NAMESPACE ./Subversion-in-BioC/OncoSimulR/.

cp OncoSimulR/DESCRIPTION ./Subversion-in-BioC/OncoSimulR/.


## should we run the tests?

if [[ $# == 1 ]]; then
    V_R=$1
    cd ~/Proyectos/OncoSimul/Subversion-in-BioC
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
cd ~/Proyectos/OncoSimul/Subversion-in-BioC/OncoSimulR

svn status --show-updates


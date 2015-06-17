#!/bin/bash
cp OncoSimulR/vignettes/OncoSimulR.Rnw ../Subversion-in-BioC/OncoSimulR/vignettes/.
cp OncoSimulR/vignettes/OncoSimulR.bib ../Subversion-in-BioC/OncoSimulR/vignettes/.
cp OncoSimulR/vignettes/gitsetinfo.sty ../Subversion-in-BioC/OncoSimulR/vignettes/.
cp OncoSimulR/vignettes/gitinfo.sty ../Subversion-in-BioC/OncoSimulR/vignettes/.
cp OncoSimulR/vignettes/gitHeadInfo.gin ../Subversion-in-BioC/OncoSimulR/vignettes/.

cp OncoSimulR/src/*.c ../Subversion-in-BioC/OncoSimulR/src/.
cp OncoSimulR/src/*.cpp ../Subversion-in-BioC/OncoSimulR/src/.
cp OncoSimulR/src/*.h ../Subversion-in-BioC/OncoSimulR/src/.

cp OncoSimulR/R/*.R ../Subversion-in-BioC/OncoSimulR/R/.

cp OncoSimulR/man/*.Rd ../Subversion-in-BioC/OncoSimulR/man/.

cp OncoSimulR/inst/NEWS ../Subversion-in-BioC/OncoSimulR/inst/.
cp OncoSimulR/inst/CITATION ../Subversion-in-BioC/OncoSimulR/inst/.
cp OncoSimulR/inst/miscell/example-binom-problems.cpp ../Subversion-in-BioC/OncoSimulR/inst/miscell/.

cp OncoSimulR/data/*.RData ../Subversion-in-BioC/OncoSimulR/data/.

cp OncoSimulR/NAMESPACE ../Subversion-in-BioC/OncoSimulR/.

cp OncoSimulR/DESCRIPTION ../Subversion-in-BioC/OncoSimulR/.


## check if things need adding

cd ~/Proyectos/OncoSimulR/Subversion-in-BioC/OncoSimulR

svn status --show-updates


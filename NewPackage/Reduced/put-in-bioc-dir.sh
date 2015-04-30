#!/bin/bash
cp ADaCGH2/vignettes/ADaCGH2.Rnw ~/Proyectos/adacgh2/R-packages/Subversion-in-BioC/ADaCGH2/vignettes/.
cp ADaCGH2/vignettes/ADaCGH2.bib ~/Proyectos/adacgh2/R-packages/Subversion-in-BioC/ADaCGH2/vignettes/.
cp ADaCGH2/inst/doc/ADaCGH2-long-examples.* ~/Proyectos/adacgh2/R-packages/Subversion-in-BioC/ADaCGH2/inst/doc/.
cp ADaCGH2/inst/doc/benchmarks.pdf ~/Proyectos/adacgh2/R-packages/Subversion-in-BioC/ADaCGH2/inst/doc/.
cp ADaCGH2/R/ADaCGH-2.R ~/Proyectos/adacgh2/R-packages/Subversion-in-BioC/ADaCGH2/R/.
cp ADaCGH2/DESCRIPTION ~/Proyectos/adacgh2/R-packages/Subversion-in-BioC/ADaCGH2/.
cp ADaCGH2/NAMESPACE ~/Proyectos/adacgh2/R-packages/Subversion-in-BioC/ADaCGH2/.
cp ADaCGH2/man/*.Rd ~/Proyectos/adacgh2/R-packages/Subversion-in-BioC/ADaCGH2/man/.
cp ADaCGH2/inst/NEWS ~/Proyectos/adacgh2/R-packages/Subversion-in-BioC/ADaCGH2/inst/.
cp ADaCGH2/inst/sources-tex-Rnw-other/* ~/Proyectos/adacgh2/R-packages/Subversion-in-BioC/ADaCGH2/inst/sources-tex-Rnw-other/.
grep -v "^%" ADaCGH2/inst/sources-tex-Rnw-other/benchmarks.tex > /home/ramon/Proyectos/adacgh2/R-packages/Subversion-in-BioC/ADaCGH2/inst/sources-tex-Rnw-other/benchmarks.tex
cp ADaCGH2/inst/doc/index.html ~/Proyectos/adacgh2/R-packages/Subversion-in-BioC/ADaCGH2/inst/doc/.
cp ADaCGH2/README.code.authors ~/Proyectos/adacgh2/R-packages/Subversion-in-BioC/ADaCGH2/.
cp ADaCGH2/src/* ~/Proyectos/adacgh2/R-packages/Subversion-in-BioC/ADaCGH2/src/.


cp ADaCGH2/inst/CITATION ~/Proyectos/adacgh2/R-packages/Subversion-in-BioC/ADaCGH2/inst/.


## check if things need adding

## data not coppied, since very rarely changed

cd ~/Proyectos/adacgh2/R-packages/Subversion-in-BioC/ADaCGH2

svn status --show-updates

## to create the directory
## tar --exclude=.svn -zcvf ADaCGH2-como-esta-en-BioC.tar.gz ADaCGH2

#!/bin/bash
## simplify my workflow.
cd ./OncoSimulR/vignettes
Rscript-devel -e 'library(rmarkdown); library(BiocStyle); render("OncoSimulR.Rmd")'
cd ../..

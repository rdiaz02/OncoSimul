#!/bin/bash
## simplify my workflow.
cd ./OncoSimulR/vignettes
Rscript -e 'library(rmarkdown); library(BiocStyle); render("OncoSimulR.Rmd")'
cd ../..

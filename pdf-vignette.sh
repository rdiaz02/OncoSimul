#!/bin/bash
## simplify my workflow.
cd ./OncoSimulR/vignettes
Rscript -e 'library(rmarkdown); library(BiocStyle); library(bookdown); render("OncoSimulR.Rmd", output_format = bookdown::pdf_document2(toc = TRUE, toc_depth = 4, keep_tex = TRUE))'
cd ../..

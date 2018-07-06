#! /bin/bash
## for quick things, faster than R CMD check, etc
R -e 'library(knitr); knit("../OncoSimulR/vignettes/OncoSimulR.Rmd")'
R -e 'library(testthat); library(gtools); library(smatr); test_dir("../OncoSimulR/tests/testthat")'
R -e 'library(testthat); library(gtools); library(smatr); test_dir("../OncoSimulR/tests/manual")'

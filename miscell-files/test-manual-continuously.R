## For running the manual tests

library(testthat)
library(OncoSimulR)
library(help = "OncoSimulR")
library(car)
library(smatr)
library(gtools)
## recall to install with --install-tests

if (Sys.info()["nodename"] == "Gallotia") {
    td <- "/home2/ramon/tmp"
} else {
    td <- "~/tmp"
}


i <- 0
while(TRUE) {
    i <- i + 1
    cat("\n\n Doing iteration ", i, "\n")
    cat("\n And this is the first random uniform number ", runif(1), "\n")
    the.seed <- .Random.seed ## tests set the seed in several places
    gc()
    test_dir("../OncoSimulR/tests/manual/")
    gc()
    .Random.seed <- the.seed
    cat("\n         And this is the second random uniform number ", runif(1), "\n")
}



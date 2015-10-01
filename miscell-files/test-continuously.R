## Recall we can just use R CMD check

library(testthat)
library(OncoSimulR)
library(tools)
## recall to install with --install-tests

i <- 0
while(TRUE) {
    i <- i + 1
    cat("\n\n Doing iteration ", i, "\n")
    cat("\n And this is a random uniform number ", runif(1), "\n")
    nn <- paste(sample(c(letters, 0:1), 12, replace = TRUE), collapse = "")
    nnn <- paste0(tempfile(pattern=""), nn)
    dir.create(nnn)
    cat(paste("Using directory ", nnn, "\n"))
    test_package("OncoSimulR") ## testthat, but if interactive, asks for plots,
                               ## and you need to load the library first
    testInstalledPackage(pkg = "OncoSimulR", outDir = nnn,
                         types = c("examples", "tests"))
    checkVignettes("OncoSimulR", workdir = "src") ## yes, tmp or cur does not do it
    cat("\n And this is the final random uniform number ", runif(1), "\n")
}



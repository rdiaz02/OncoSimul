## Recall we can just use R CMD check

library(testthat)
library(OncoSimulR)
library(tools)
## recall to install with --install-tests

i <- 0
while(TRUE) {
    i <- i + 1
    cat("\n\n Doing iteration ", i, "\n")
    cat("\n And this is the first random uniform number ", runif(1), "\n")
    nn <- paste(sample(c(letters, 0:1), 12, replace = TRUE), collapse = "")
    nnn <- paste0(tempfile(pattern=""), nn)
    dir.create(nnn)
    cat(paste("Using directory ", nnn, "\n"))
    test_package("OncoSimulR") ## testthat, but if interactive, asks for plots,
    ## and you need to load the library first
    the.seed <- .Random.seed ## examples and vignette set the seed in several
    ## places
    cat("\n         And this is the second random uniform number ", runif(1), "\n")

    testInstalledPackage(pkg = "OncoSimulR", outDir = nnn,
                         types = c("examples"))
    cat("\n                    And this is the third random uniform number ", runif(1), "\n")
    .Random.seed <- the.seed
    testInstalledPackage(pkg = "OncoSimulR", outDir = nnn,
                         types = c("tests"))
    cat("\n                          And this is the fourth random uniform number ", runif(1), "\n")
    ## .Random.seed <- the.seed  ## note the seed does not advance here at all. Did above.
    ## checkVignettes("OncoSimulR", workdir = "src") ## yes, tmp or cur does not do it
    cat("\n                              And this is the FINAL random uniform number ", runif(1), "\n")
    .Random.seed <- the.seed
}



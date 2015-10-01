## Recall we can just use R CMD check

library(testthat)
library(OncoSimulR)
library(tools)
## recall to install with --install-tests

while(TRUE) {
    nn <- paste(sample(c(letters, 0:1), 12, replace = TRUE), collapse = "")
    nnn <- paste0("~/tmp/", nn)
    dir.create(nnn)
    cat(paste("Using directory ", nnn, "\n"))
    test_package("OncoSimulR") ## testthat, but if interactive, asks for plots,
                               ## and you need to load the library first
    testInstalledPackage(pkg = "OncoSimulR", outDir = nnn,
                         types = c("examples", "tests"))
    checkVignettes("OncoSimulR", workdir = "src") ## yes, tmp or cur does not do it
}



## Recall we can just use R CMD check

library(testthat)
library(OncoSimulR)
library(tools)
## library(devtools)
library(knitr)

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
    nn <- paste(sample(c(letters, 0:1), 15, replace = TRUE), collapse = "")
    nnn <- paste0(tempfile(pattern="Rtmp", tmpdir = td), nn)
    dir.create(nnn)
    cat(paste("Using directory ", nnn, "\n"))

    cat(paste("\n starting test_package", date()))
    test_package("OncoSimulR") ## testthat, but if interactive, asks for plots,
    ## and you need to load the library first
    cat(paste("\n  done test_package", date(), "\n"))
    
    cat("\n         And this is the second random uniform number ", runif(1), "\n")

    cat(paste("\n starting testInstalledPackage, examples", date()))
    testInstalledPackage(pkg = "OncoSimulR", outDir = nnn,
                         types = c("examples"))
    cat(paste("\n  done testInstalledPackage, examples", date(), "\n"))

    cat("\n                    And this is the third random uniform number ", runif(1), "\n")
    cat(paste("\n starting testInstalledPackage, tests", date()))
    testInstalledPackage(pkg = "OncoSimulR", outDir = nnn,
                         types = c("tests"))
    cat(paste("\n  done testInstalledPackage, tests", date(), "\n"))

    cat("\n                          And this is the fourth random uniform number ", runif(1), "\n")

    ## Vignette via knitr.
    nf <- tempfile()
    cat(paste("\n     knit output to ", nf, "\n"))
    
    knit("../OncoSimulR/vignettes/OncoSimulR.Rnw", output = nf)

    ## if you want to tex the file, use knit2pdf and change output name to
    ## have tex extension, etc.

    ## This if using the source But either I do not understand or this
    ## does not work, as using workdir = "tmp" fails. And what will happen
    ## with multiple parallel jobs?
    ## checkVignettes(dir = "../OncoSimulR", workdir = "src") ## yes, tmp or
    ##                                                        ## cur does not
    ##                                                        ## do it

    
    ## devtools, checks vignette and also catches other errors. ALWAYS
    ## commit changes before as devtools can change files!!!  But the
    ## changing of files is way too dangerous for a check.
    ## devtools::check(pkg = "../OncoSimulR", document = FALSE)

    ## cat(paste("\n starting testInstalledPackage, vignette", date()))
    ## This does not work. 
    ## testInstalledPackage(pkg = "OncoSimulR", outDir = nnn, types =
    ## c("vignettes")) cat(paste("\n done testInstalledPackage, vignette",
    ## date(), "\n"))
    ## checkVignettes("OncoSimulR", workdir = "src") ## yes, tmp or cur does not do it
    cat("\n                              And this is the FINAL random uniform number ", runif(1), "\n")
    ## set.seed(NULL) ## Not needed as files clean up after themselves
}



## Recall we can just use R CMD check

## FIXME I do not want to keep track of all the output files but cannot
## simply write to /dev/null. Check https://github.com/xrgtn/nullfs or
## write in /tmp and remove. Have /tmp mounted in tmpfs, in RAM. Myabe
## this is the simplest.

library(testthat)
library(OncoSimulR)
library(tools)
## library(devtools)
library(knitr)

## recall to install with --install-tests

## if (Sys.info()["nodename"] == "Gallotia") {
##     td <- "/home2/ramon/tmp"
## } else {
##     td <- "~/tmp"
## }

## If /tmp is tmpfs in RAM much faster. And we clean up. So use that.
td <- "/tmp"

iter <- 0
while(TRUE) {
    iter <- iter + 1
    cat("\n\n Doing iteration ", iter, "\n")
    cat("\n And this is the first random uniform number ", runif(1), "\n")
    nn <- paste(sample(c(letters, 0:1), 15, replace = TRUE), collapse = "")
    nnn <- paste0(tempfile(pattern="Rtmp", tmpdir = td), nn)
    dir.create(nnn)
    cat(paste("Using directory ", nnn, "\n"))

    cat(paste("\n starting test_package", date()))
    test_package("OncoSimulR") ## testthat, but if interactive, asks for plots,
    ## and you need to load the library first
    cat(paste("\n  done test_package", date(), "\n"))
    gc(); gc(); gc(); gc()
    
    cat("\n         And this is the second random uniform number ", runif(1), "\n")

    cat(paste("\n starting testInstalledPackage, examples", date()))
    testInstalledPackage(pkg = "OncoSimulR", outDir = nnn,
                         types = c("examples"))
    cat(paste("\n  done testInstalledPackage, examples", date(), "\n"))
    gc(); gc(); gc(); gc()
    cat("\n                    And this is the third random uniform number ", runif(1), "\n")
    cat(paste("\n starting testInstalledPackage, tests", date()))
    testInstalledPackage(pkg = "OncoSimulR", outDir = nnn,
                         types = c("tests"))
    cat(paste("\n  done testInstalledPackage, tests", date(), "\n"))
    gc(); gc(); gc(); gc()
    cat("\n                          And this is the fourth random uniform number ", runif(1), "\n")

    ## Vignette via knitr.
    nf <- tempfile(tmpdir = nnn)
    cat(paste("\n     knit output to ", nf, "\n"))
    
    knit("../OncoSimulR/vignettes/OncoSimulR.Rmd", output = nf)
    unlink(nf, recursive = TRUE)
    gc(); gc(); gc(); gc()
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
    ## rm the tmp
    cat("\n  Cleaning up tmp files \n")
    unlink(nnn, recursive = TRUE)
}



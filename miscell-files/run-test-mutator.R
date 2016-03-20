library(OncoSimulR)
library(testthat)

while(TRUE) {
    RNGkind("Mersenne-Twister")
    cat(paste("\n a runif is ", runif(1), "\n"))
    the.seed <- .Random.seed
    try(source("../OncoSimulR/tests/testthat/test.mutator.R", echo = TRUE,
               max.deparse.length = 9999999999999))
    .Random.seed <- the.seed
    cat(paste("\n    Two runifs are ", runif(2), "\n"))
}

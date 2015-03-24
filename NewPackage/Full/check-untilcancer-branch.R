## We run the following script with two versions, 99.1.4, which is from
## before the "untilcancer" in cpp, and 99.1.5, from after. In CPP we set
## the DEBUGV, so very verbose. And compare. They should be identical.


## run as
## R-patched --vanilla --args 5 < check-untilcancer-branch.R > check-untilcancer-branch-untilcancer.Rout

args <- commandArgs(TRUE)
if(args  == "4") {
    ## this is where I install the previous library
    library(OncoSimulR, lib.loc = "/home/ramon/tmp")
} else {
    library(OncoSimulR)
}

packageVersion("OncoSimulR")

## I set onlyCancer = FALSE, because repeated iterations, in the first
## version, go back to R, which then generates a new seed for C++.

data(examplePosets)
p701 <- examplePosets[["p701"]]
set.seed(1)
b1 <- oncoSimulIndiv(p701, verbosity = 3,
                     initSize = 500,
                     detectionSize = 700,
                     onlyCancer = FALSE)


set.seed(25)
m1 <- oncoSimulIndiv(p701,
                     model = "McFL",
                     mu = 5e-7,
                     initSize = 4000,
                     detectionDrivers = 1,
                     sampleEvery = 0.025, # so it is faster
                     finalTime = 15000,
                     endTimeEvery = 5 * 0.025,
                     keepEvery = 5, verbosity = 3,
                     onlyCancer = FALSE)

set.seed(67)
m2 <- oncoSimulIndiv(p701,
                     model = "McFL",
                     mu = 5e-5,
                     initSize = 8000,
                     detectionDrivers = 1,
                     sampleEvery = 0.025, # so it is faster
                     finalTime = 15000,
                     endTimeEvery = -9,
                     keepEvery = 5, verbosity = 3,
                     onlyCancer = FALSE)

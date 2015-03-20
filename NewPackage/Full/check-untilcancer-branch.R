## We run the following script with two versions, 99.1.4, which is from
## before the "untilcancer" in cpp, and 99.1.5, from after. In CPP we set
## the DEBUGV, so very verbose. And compare. They should be identical.

library(OncoSimulR)
packageVersion("OncoSimulR")

data(examplePosets)
p701 <- examplePosets[["p701"]]
set.seed(1)
b1 <- oncoSimulIndiv(p701, verbosity = 3,
                     initSize = 500,
                     detectionSize = 2000,
                     onlyCancer = FALSE)


set.seed(1)
m1 <- oncoSimulIndiv(p701,
                     model = "McFL",
                     mu = 5e-7,
                     initSize = 4000,
                     detectionDrivers = 1,
                     sampleEvery = 0.025,
                     finalTime = 15000,
                     keepEvery = 5, verbosity = 3,
                     onlyCancer = FALSE)

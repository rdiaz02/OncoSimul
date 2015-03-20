## R-patched CMD build --no-build-vignettes OncoSimulR; R-patched CMD INSTALL OncoSimulR_99.1.4.tar.gz

library(OncoSimulR)
data(examplePosets)
p1101 <- examplePosets[["p1101"]]
b1 <- oncoSimulIndiv(p1101, keepEvery = 15, verbosity = 0)


m1 <- oncoSimulIndiv(p1101,
                     model = "McFL",
                     mu = 5e-7,
                     initSize = 4000,
                     sampleEvery = 0.025,
                     finalTime = 15000,
                     keepEvery = 5,
                     verbosity = 3)

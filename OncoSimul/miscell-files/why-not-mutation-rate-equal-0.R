We could if we wanted to, but it makes no sense. And we would need to add
lots of code to handle at least two special cases:

- That we can end up with mutable positions (positions that are not yet
  mutated) but have a mutation rate of 0. In general, that is not allowd
  and triggers and exception.

- That we would see the gene with mu = 0 eventually get mutated as soon as
  we start either: a) reaching the dummy mutation scenario; b) just using
  a discrete_distribution in C++ where there is only one item; evne if the
  weight is 0, it gets, of course, chosen. Use ex-dics-dist.cpp and play
  with it. There is no way around that.







## This breaks with Unrecoverable exception: mutation = 0 with numMutable
## != 0?. Aborting.
## And that makes sense, of course.
## set.seed(26)
## muvar2 <- c("U" = 0, "D" = 1e-3)
## ni2 <- rep(0.01, length(muvar2))
## names(ni2) <- names(muvar2)
## fe1 <- allFitnessEffects(noIntGenes = ni2)
## no <- 1e4
## reps <- 100
## bb <- oncoSimulIndiv(fe1, mu = muvar2, onlyCancer = FALSE,
##                      initSize = no,
##                      model = "McFL",
##                      finalTime = 200,
##                      seed = NULL
##                      )
## bb


## Shows getting the gene with mutation = 0
## set.seed(26)
## muvar2 <- c("U" = 0, "D" = 1e-3, "e" = 1e-2)
## ni2 <- rep(0.01, length(muvar2))
## names(ni2) <- names(muvar2)
## fe1 <- allFitnessEffects(noIntGenes = ni2)
## no <- 1e4
## reps <- 100
## bb <- oncoSimulIndiv(fe1, mu = muvar2, onlyCancer = FALSE,
##                      initSize = no,
##                      model = "McFL",
##                      finalTime = 200,
##                      verbosity = 3,
##                      seed = NULL
##                      )
## bb

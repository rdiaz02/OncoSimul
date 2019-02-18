## this can do
sd <- 0.1 ## fitness effect of drivers
sm <- 0 ## fitness effect of mutator
nd <- 20 ## number of drivers
nm <- 5  ## number of mutators
mut <- 50 ## mutator effect

## Blows up because the mutation rate gets huge

fitnessGenesVector <- c(rep(sd, nd), rep(sm, nm))
names(fitnessGenesVector) <- 1:(nd + nm)
mutatorGenesVector <- rep(mut, nm)
names(mutatorGenesVector) <- (nd + 1):(nd + nm)

ft <- allFitnessEffects(noIntGenes = fitnessGenesVector,
                        drvNames = 1:nd)
mt <- allMutatorEffects(noIntGenes = mutatorGenesVector)
set.seed(2)
RNGkind("L'Ecuyer-CMRG")
st <- oncoSimulPop(4, ft, muEF = mt,
                   detectionDrivers = 4,
                   finalTime = NA,
                   detectionSize = NA,
                   detectionProb = NA,
                   onlyCancer = TRUE,
                   mc.cores = 2, ## adapt to your hardware
                   seed = NULL)

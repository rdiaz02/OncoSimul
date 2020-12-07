set.seed(1)
pops <- 200
ft <- 5e-3
lni <- 7
no <- 5e5
ni <- c(0, 0, 0, rep(0, lni))
## scramble around names
names(ni) <- c("hereisoneagene",
               "oreoisasabgene",
               "nnhsisthecgene",
               replicate(lni,
                         paste(sample(letters, 12), collapse = "")))
ni <- ni[order(names(ni))]
fe <- allFitnessEffects(noIntGenes = ni)
mutator1 <- rep(1, lni + 3)
pg1 <- seq(from = 1e-9, to = 1e-6, length.out = lni + 3) ## max should not be
## huge here as mutator
## is 34. Can get beyond
## 1
names(mutator1) <- sample(names(ni))
names(pg1) <- sample(names(ni))
mutator1["oreoisasabgene"] <- 100
m1 <- allMutatorEffects(noIntGenes = mutator1)
set.seed(1)
oncoSimulIndiv(fe,  detectionProb = NA,
               mu = pg1,
               muEF = m1,
               finalTime = ft,
               mutationPropGrowth = FALSE,
               initSize = no,
               initMutant ="oreoisasabgene",
               model = "McFL",
               sampleEvery = 0.01, 
               detectionSize = 1e9,
               detectionDrivers = 9999,
               onlyCancer = FALSE, seed = NULL, verbosity = 5)

## This crashes
set.seed(1) ## for reproducibility
## 17 genes, 7 with no direct fitness effects
ni <- c(rep(0, 7), runif(10, min = -0.01, max = 0.1))
names(ni) <- c("a1", "a2", "b1", "b2", "b3", "c1", "c2",
               paste0("g", 1:10))
fe31 <- allFitnessEffects(noIntGenes = ni,
                         drvNames = "g1")
drvNames = paste0("g", 1:10)
set.seed(1) ## so that it is easy to reproduce
mue11 <- oncoSimulIndiv(fe31, 
                       mu = 1e-6,
                       initSize = 1e5,
                       model = "McFL",
                       detectionSize = 5e6,
                       finalTime = 2,
                       onlyCancer = FALSE)



## This does not
set.seed(1) ## for reproducibility
## 17 genes, 7 with no direct fitness effects
ni <- c(rep(0, 7), runif(10, min = -0.01, max = 0.1))
names(ni) <- c("a1", "a2", "b1", "b2", "b3", "c1", "c2",
               paste0("g", 1:10))
fe31 <- allFitnessEffects(noIntGenes = ni,
                         drvNames = "a2")
drvNames = paste0("g", 1:10)
set.seed(1) ## so that it is easy to reproduce
mue11 <- oncoSimulIndiv(fe31, 
                       mu = 1e-6,
                       initSize = 1e5,
                       model = "McFL",
                       detectionSize = 5e6,
                       finalTime = 2,
                       onlyCancer = FALSE)

## Problem is having a noIntGene as a driver




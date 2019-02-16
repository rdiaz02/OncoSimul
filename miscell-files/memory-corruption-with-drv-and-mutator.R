## This crashes
set.seed(1) ## for reproducibility
## 17 genes, 7 with no direct fitness effects
ni <- c(rep(0, 7), runif(10, min = -0.01, max = 0.1))
names(ni) <- c("a1", "a2", "b1", "b2", "b3", "c1", "c2",
               paste0("g", 1:10))
fe3 <- allFitnessEffects(noIntGenes = ni,
                         drvNames = paste0("g", 1:10))
fm3 <- allMutatorEffects(epistasis = c("A" = 5,
                                       "B" = 10,
                                       "C" = 3,
                                       "A:C" = 70),
                         geneToModule = c("A" = "a1, a2",
                                          "B" = "b1, b2, b3",
                                          "C" = "c1, c2"))
evalGenotypeFitAndMut("a1, a2", fe3, fm3, verbose = TRUE)
evalGenotypeFitAndMut("a1, b3", fe3, fm3, verbose = TRUE)
## These only affect fitness: the mutator multiplier is 1
evalGenotypeFitAndMut("g1", fe3, fm3, verbose = TRUE)                      
evalGenotypeFitAndMut("g3, g9", fe3, fm3, verbose = TRUE)
set.seed(1) ## so that it is easy to reproduce
mue1 <- oncoSimulIndiv(fe3, muEF = fm3, 
                       mu = 1e-6,
                       initSize = 1e5,
                       model = "McFL",
                       detectionSize = 5e6,
                       finalTime = 1,
                       onlyCancer = FALSE)




## The next is OK
set.seed(1) ## for reproducibility
## 17 genes, 7 with no direct fitness effects
ni <- c(rep(0, 7), runif(10, min = -0.01, max = 0.1))
names(ni) <- c("a1", "a2", "b1", "b2", "b3", "c1", "c2",
               paste0("g", 1:10))
fe3 <- allFitnessEffects(noIntGenes = ni,
                         drvNames = c("a1", "a2"))
fm3 <- allMutatorEffects(epistasis = c("A" = 5,
                                       "B" = 10,
                                       "C" = 3,
                                       "A:C" = 70),
                         geneToModule = c("A" = "a1, a2",
                                          "B" = "b1, b2, b3",
                                          "C" = "c1, c2"))
evalGenotypeFitAndMut("a1, a2", fe3, fm3, verbose = TRUE)
evalGenotypeFitAndMut("a1, b3", fe3, fm3, verbose = TRUE)
## These only affect fitness: the mutator multiplier is 1
evalGenotypeFitAndMut("g1", fe3, fm3, verbose = TRUE)                      
evalGenotypeFitAndMut("g3, g9", fe3, fm3, verbose = TRUE)
set.seed(1) ## so that it is easy to reproduce
mue1 <- oncoSimulIndiv(fe3, muEF = fm3, 
                       mu = 1e-6,
                       initSize = 1e5,
                       model = "McFL",
                       detectionSize = 5e6,
                       finalTime = 1,
                       onlyCancer = FALSE)



## This does not. Why??!!!
set.seed(1) ## for reproducibility
## 17 genes, 7 with no direct fitness effects
ni <- c(rep(0, 7), runif(10, min = -0.01, max = 0.1))
names(ni) <- c("a1", "a2", "b1", "b2", "b3", "c1", "c2",
               paste0("g", 1:10))
fe31 <- allFitnessEffects(noIntGenes = ni,
                         drvNames = "a2")
set.seed(1) ## so that it is easy to reproduce
mue11 <- oncoSimulIndiv(fe31, 
                       mu = 1e-6,
                       initSize = 1e5,
                       model = "McFL",
                       detectionSize = 5e6,
                       finalTime = 5,
                       onlyCancer = FALSE)

## This does not. Why??!!!
set.seed(1) ## for reproducibility
## 17 genes, 7 with no direct fitness effects
ni <- c(rep(0, 7), runif(10, min = -0.01, max = 0.1))
names(ni) <- c("a1", "a2", "b1", "b2", "b3", "c1", "c2",
               paste0("g", 1:10))
fe31 <- allFitnessEffects(noIntGenes = ni,
                         drvNames = "b3")
set.seed(1) ## so that it is easy to reproduce
mue11 <- oncoSimulIndiv(fe31, 
                       mu = 1e-6,
                       initSize = 1e5,
                       model = "McFL",
                       detectionSize = 5e6,
                       finalTime = 15,
                       onlyCancer = FALSE)
mue11


## This does Why??!!!
set.seed(1) ## for reproducibility
## 17 genes, 7 with no direct fitness effects
ni <- c(rep(0, 7), runif(10, min = -0.01, max = 0.1))
names(ni) <- c("a1", "a2", "b1", "b2", "b3", "c1", "c2",
               paste0("g", 1:10))
fe31 <- allFitnessEffects(noIntGenes = ni,
                         drvNames = "g1")
set.seed(1) ## so that it is easy to reproduce
mue11 <- oncoSimulIndiv(fe31, 
                       mu = 1e-6,
                       initSize = 1e5,
                       model = "McFL",
                       detectionSize = 5e6,
                       finalTime = 15,
                       onlyCancer = FALSE)
mue11



###

## This does not. Why??!!!
set.seed(1) ## for reproducibility
## 17 genes, 7 with no direct fitness effects
ni <- c(rep(0, 7), runif(10, min = -0.01, max = 0.1))
names(ni) <- c("g1", "a2", "b1", "b2", "b3", "c1", "c2",
               paste0("g", 2:11))
fe31 <- allFitnessEffects(noIntGenes = ni,
                         drvNames = "g1")
set.seed(1) ## so that it is easy to reproduce
mue11 <- oncoSimulIndiv(fe31, 
                       mu = 1e-6,
                       initSize = 1e5,
                       model = "McFL",
                       detectionSize = 5e6,
                       finalTime = 15,
                       onlyCancer = FALSE)
mue11


## This does not. Why??!!!
set.seed(1) ## for reproducibility
## 17 genes, 7 with no direct fitness effects
ni <- c(rep(0.1, 7), runif(10, min = -0.01, max = 0.1))
names(ni) <- c("g1", "a2", "b1", "b2", "b3", "c1", "c2",
               paste0("g", 2:11))
fe31 <- allFitnessEffects(noIntGenes = ni,
                         drvNames = "c1")
set.seed(1) ## so that it is easy to reproduce
mue11 <- oncoSimulIndiv(fe31, 
                       mu = 1e-6,
                       initSize = 1e5,
                       model = "McFL",
                       detectionSize = 5e6,
                       finalTime = 25,
                       onlyCancer = FALSE)
mue11


## This does. Why??!!!
set.seed(1) ## for reproducibility
## 17 genes, 7 with no direct fitness effects
ni <- c(rep(0.1, 7), runif(10, min = -0.01, max = 0.1))
names(ni) <- c("g1", "a2", "b1", "b2", "b3", "c1", "c2",
               paste0("g", 2:11))
fe31 <- allFitnessEffects(noIntGenes = ni,
                         drvNames = "c2")
set.seed(1) ## so that it is easy to reproduce
mue11 <- oncoSimulIndiv(fe31, 
                       mu = 1e-6,
                       initSize = 1e5,
                       model = "McFL",
                       detectionSize = 5e6,
                       finalTime = 15,
                       onlyCancer = FALSE)
mue11




## This does not. Why??!!!
set.seed(1) ## for reproducibility
## 17 genes, 7 with no direct fitness effects
ni <- runif(10, min = -0.01, max = 0.1)
names(ni) <- paste0("g", 2:11)
fe31 <- allFitnessEffects(noIntGenes = ni,
                          drvNames = "g2")
set.seed(1) ## so that it is easy to reproduce
mue11 <- oncoSimulIndiv(fe31, 
                       mu = 1e-6,
                       initSize = 1e5,
                       model = "McFL",
                       detectionSize = 5e6,
                       finalTime = 25,
                       onlyCancer = FALSE)
mue11



## This does. Why??!!!
set.seed(1) ## for reproducibility
## 17 genes, 7 with no direct fitness effects
ni <- runif(10, min = -0.01, max = 0.1)
names(ni) <- paste0("g", 2:11)
fe31 <- allFitnessEffects(noIntGenes = ni,
                          drvNames = "g7") ## but g7 does not
set.seed(1) ## so that it is easy to reproduce
mue11 <- oncoSimulIndiv(fe31, 
                       mu = 1e-6,
                       initSize = 1e5,
                       model = "McFL",
                       detectionSize = 5e6,
                       finalTime = 25,
                       onlyCancer = FALSE)
mue11


## This does . Why??!!!
set.seed(1) ## for reproducibility
## 17 genes, 7 with no direct fitness effects
ni <- runif(10, min = -0.01, max = 0.1)
names(ni) <- paste0("h", 11:20)
fe31 <- allFitnessEffects(noIntGenes = ni,
                          drvNames = "h19") ## h11 to 17 OK
set.seed(1) ## so that it is easy to reproduce
mue11 <- oncoSimulIndiv(fe31, 
                       mu = 1e-6,
                       initSize = 1e5,
                       model = "McFL",
                       detectionSize = 5e6,
                       finalTime = 25,
                       onlyCancer = FALSE)
mue11



## Problem is having a noIntGene as a driver




## This crashes
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
                       finalTime = 2,
                       onlyCancer = FALSE)


## But this works around the issue in terms of what we really want to do
set.seed(1) ## for reproducibility
## 17 genes, 7 with no direct fitness effects
ni <- c(rep(0, 7), runif(10, min = -0.01, max = 0.1))
names(ni) <- c("a1", "a2", "b1", "b2", "b3", "c1", "c2",
               paste0("g", 1:10))
## For plotting purposes, recognize those genes with s>0
sp <- which(ni > 0)
fe31 <- allFitnessEffects(noIntGenes = ni[-sp],
                          epistasis = ni[sp],
                          drvNames = names(ni)[sp])
set.seed(1) ## so that it is easy to reproduce
mue11 <- oncoSimulIndiv(fe31, 
                       mu = 1e-6,
                       initSize = 1e5,
                       model = "McFL",
                       detectionSize = 5e6,
                       finalTime = 2,
                       onlyCancer = FALSE)
## and then, check that the above fails, and this works,
## and epistasis is correct if genes with no interactions.


## This does not
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
                       finalTime = 2,
                       onlyCancer = FALSE)

## Problem is having a noIntGene as a driver




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
                          drvNames = "g7") ## g9 crashes, but g7 does not
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



## This does. Why??!!!
set.seed(1) ## for reproducibility
## 17 genes, 7 with no direct fitness effects
ni <- runif(10, min = -0.01, max = 0.1)
names(ni) <- paste0("g", 2:11)
fe31 <- allFitnessEffects(noIntGenes = ni,
                          drvNames = "g7") ## g9 crashes, but g7 does not
set.seed(1) ## so that it is easy to reproduce
mue11 <- oncoSimulIndiv(fe31, 
                       mu = 1e-6,
                       initSize = 1e5,
                       model = "McFL",
                       detectionSize = 5e6,
                       finalTime = 2.8,
                       onlyCancer = FALSE)
mue11


## This does not. Why??!!!
set.seed(1) ## for reproducibility
## 17 genes, 7 with no direct fitness effects
ni <- runif(10, min = -0.01, max = 0.1)
names(ni) <- paste0("g", 2:11)
fe31 <- allFitnessEffects(noIntGenes = ni,
                          drvNames = "g3") ## g9 crashes, but g7 does not
set.seed(1) ## so that it is easy to reproduce
mue11 <- oncoSimulIndiv(fe31, 
                       mu = 1e-6,
                       initSize = 1e5,
                       model = "McFL",
                       detectionSize = 5e6,
                       finalTime = 2.8,
                       onlyCancer = FALSE)
mue11


## This does. Why??!!!
set.seed(1) ## for reproducibility
## 17 genes, 7 with no direct fitness effects
ni <- runif(10, min = -0.01, max = 0.1)
names(ni) <- paste0("g", 2:11)
fe31 <- allFitnessEffects(noIntGenes = ni,
                          drvNames = c("g3", "g8")) 
set.seed(1) ## so that it is easy to reproduce
mue11 <- oncoSimulIndiv(fe31, 
                       mu = 1e-6,
                       initSize = 1e5,
                       model = "McFL",
                       detectionSize = 5e6,
                       finalTime = 2.8,
                       onlyCancer = FALSE)
mue11


## This does not. Why??!!!
set.seed(1) ## for reproducibility
## 17 genes, 7 with no direct fitness effects
ni <- runif(10, min = -0.01, max = 0.1)
names(ni) <- paste0("g", 2:11)
fe31 <- allFitnessEffects(noIntGenes = ni,
                          drvNames = names(ni)) ## g9 crashes, but g7 does not
set.seed(1) ## so that it is easy to reproduce
mue11 <- oncoSimulIndiv(fe31, 
                       mu = 1e-6,
                       initSize = 1e5,
                       model = "McFL",
                       detectionSize = 5e6,
                       finalTime = 2.7,
                       onlyCancer = FALSE)
mue11


## This does not. Why??!!!
set.seed(1) ## for reproducibility
## 17 genes, 7 with no direct fitness effects
ni <- runif(10, min = -0.01, max = 0.1)
names(ni) <- paste0("g", 2:11)
fe31 <- allFitnessEffects(noIntGenes = ni,
                          drvNames = names(ni)) ## g9 crashes, but g7 does not
set.seed(1) ## so that it is easy to reproduce
mue11 <- oncoSimulIndiv(fe31, 
                       mu = 1e-6,
                       initSize = 1e5,
                       model = "McFL",
                       detectionSize = 5e6,
                       finalTime = 2.8,
                       onlyCancer = FALSE)
mue11


## This does . Why??!!!
set.seed(1) ## for reproducibility
## 17 genes, 7 with no direct fitness effects
ni <- runif(10, min = -0.01, max = 0.1)
names(ni) <- paste0("g", 2:11)
fe31 <- allFitnessEffects(noIntGenes = ni,
                          drvNames = names(ni)[8]) ## g9 crashes, but g7 does not
set.seed(1) ## so that it is easy to reproduce
mue11 <- oncoSimulIndiv(fe31, 
                       mu = 1e-6,
                       initSize = 1e5,
                       model = "McFL",
                       detectionSize = 5e6,
                       finalTime = 2.8,
                       verbosity = 7,
                       onlyCancer = FALSE)
mue11





## This does . Why??!!!
set.seed(1) ## for reproducibility
## 17 genes, 7 with no direct fitness effects
ni <- runif(10, min = -0.01, max = 0.1)
names(ni) <- paste0("g", 2:11)
fe31 <- allFitnessEffects(noIntGenes = ni,
                          drvNames = names(ni)[8]) ## g9 crashes, but g7 does not
set.seed(1) ## so that it is easy to reproduce
mue11 <- oncoSimulIndiv(fe31, 
                       mu = 1e-6,
                       initSize = 1e5,
                       model = "McFL",
                       detectionSize = 5e6,
                       finalTime = 2.48,
                       onlyCancer = FALSE)
mue11



## This does . Why??!!!
set.seed(1) ## for reproducibility
## 17 genes, 7 with no direct fitness effects
ni <- runif(10, min = -0.01, max = 0.1)
names(ni) <- paste0("g", 2:11)
fe31 <- allFitnessEffects(noIntGenes = ni,
                          drvNames = names(ni)[8]) ## g9 crashes, but g7 does not
set.seed(1) ## so that it is easy to reproduce
mue11 <- oncoSimulIndiv(fe31, 
                       mu = 1e-6,
                       initSize = 1e5,
                       model = "McFL",
                       detectionSize = 5e6,
                       finalTime = 2.5,
                       ## verbosity = 4, ## a 4  prevents the bug!!
                       onlyCancer = FALSE)
mue11


.Call("OncoSimulR_nr_BNB_Algo5", PACKAGE = "OncoSimulR", 
        rFE, mu, death, initSize, sampleEvery, detectionSize, 
        finalTime, initSp, initIt, seed, verbosity, speciesFS, 
        ratioForce, typeFitness_, maxram, mutationPropGrowth, 
        initMutant_, maxWallTime, keepEvery, alpha, K, detectionDrivers, 
        onlyCancer, errorHitWallTime, maxNumTries, errorHitMaxTries, 
        minDetectDrvCloneSz, extraTime, keepPhylog)











set.seed(456)
nd <- 70  
np <- 5000 
s <- 0.1  
sp <- 1e-3 
spp <- -sp/(1 + sp)
mcf1 <- allFitnessEffects(noIntGenes = c(rep(s, nd), rep(spp, np)),
                          drvNames = seq.int(nd))
mcf1s <-  oncoSimulIndiv(mcf1,
                         model = "McFL", 
                         mu = 1e-7,
                         detectionSize = 1e8, 
                         detectionDrivers = 100,
                         sampleEvery = 0.02,
                         keepEvery = 2,
                         initSize = 2000,
                         finalTime = 1000,
                         onlyCancer = FALSE)


## this is OK too
set.seed(456)
nd <- 70  
np <- 5000 
s <- 0.1  
sp <- 1e-3 
spp <- -sp/(1 + sp)
ni <- c(rep(s, nd), rep(spp, np))
names(ni) <- paste0("ng", seq.int(length(ni)))
mcf1 <- allFitnessEffects(noIntGenes = ni,
                          drvNames = names(ni))
mcf1s <-  oncoSimulIndiv(mcf1,
                         model = "McFL", 
                         mu = 1e-7,
                         detectionSize = 1e8, 
                         detectionDrivers = 100,
                         sampleEvery = 0.02,
                         keepEvery = 2,
                         initSize = 2000,
                         finalTime = 1000,
                         onlyCancer = FALSE)




## This breaks 
set.seed(1) ## for reproducibility
## 17 genes, 7 with no direct fitness effects
ni <- runif(10, min = -0.01, max = 0.1)
names(ni) <- paste0("hmuv", 11:20)
fe31 <- allFitnessEffects(noIntGenes = ni,
                          drvNames = c("hmuv19")) ## h11 to 17 OK
set.seed(1) ## so that it is easy to reproduce
mue11 <- oncoSimulIndiv(fe31, 
                       mu = 1e-6,
                       initSize = 1e5,
                       model = "McFL",
                       detectionSize = 5e6,
                       finalTime = 25,
                       onlyCancer = FALSE)
mue11


## This crashes
set.seed(1) ## for reproducibility
## 17 genes, 7 with no direct fitness effects
ni <- c(rep(0, 7), runif(10, min = -0.01, max = 0.1))
names(ni) <- c("a1", "a2", "b1", "b2", "b3", "c1", "c2",
               paste0("g", 1:10))
fe31 <- allFitnessEffects(noIntGenes = ni,
                         drvNames = c("g1"))
set.seed(1) ## so that it is easy to reproduce
mue11 <- oncoSimulIndiv(fe31, 
                       mu = 1e-6,
                       initSize = 1e5,
                       model = "McFL",
                       detectionSize = 5e6,
                       finalTime = 2,
                       onlyCancer = FALSE)







## This does not. Why??!!!
rm(list = ls())
set.seed(1) ## for reproducibility
ng <- 10
## 17 genes, 7 with no direct fitness effects
ni <- runif(ng, min = -0.01, max = 0.1)
names(ni) <- paste0("g", 1:ng)
fe31 <- allFitnessEffects(noIntGenes = ni,
                          drvNames = "g9") ## g9 crashes, but g7 does not
set.seed(1) ## so that it is easy to reproduce
mue11 <- oncoSimulIndiv(fe31, 
                       mu = 1e-6,
                       initSize = 1e7,
                       model = "McFL",
                       detectionSize = 5e6,
                       finalTime = 3,
                       onlyCancer = FALSE)
mue11








## This crashes
set.seed(1) ## for reproducibility
## 17 genes, 7 with no direct fitness effects
ni <- c(rep(0, 7), runif(10, min = -0.01, max = 0.1))
names(ni) <- c("a1", "a2", "b1", "b2", "b3", "c1", "c2",
               paste0("g", 1:10))
fe31 <- allFitnessEffects(noIntGenes = ni,
                         drvNames = c("g1"))
set.seed(1) ## so that it is easy to reproduce
mue11 <- oncoSimulIndiv(fe31, 
                       mu = 1e-6,
                       initSize = 1e6,
                       model = "McFL",
                       detectionSize = 5e6,
                       finalTime = 2,
                       onlyCancer = FALSE)










#while(TRUE){

## This does. Why??!!!
rm(list = ls())
set.seed(1) ## for reproducibility
## 17 genes, 7 with no direct fitness effects
ni <- runif(10, min = -0.01, max = 0.1)
names(ni) <- paste0("g", 2:11)
fe31 <- allFitnessEffects(noIntGenes = ni,
                          drvNames = "g9") ## g9 and sometimes g8 crashes, but g7 does not
set.seed(1) ## so that it is easy to reproduce
mue11 <- oncoSimulIndiv(fe31, 
                       mu = 1e-6,
                       initSize = 1e5,
                       model = "McFL",
                       detectionSize = 5e6,
                       finalTime = 5, # 2.43 breaks it
                       ## verbosity = 5,
                       onlyCancer = FALSE)
mue11
##}



## Moreover, with the above, we should only get maxDriversLast of 1 when
## drivers are 3 or 8 or 10; and never a 2; and MaxNumDrivers same here:
## only 1. Etc. Prepare tests that check this.


## Generate data where there are two genotypes, and consider them
## sometimes drivers or not, and make sure MaxDriversLast match.












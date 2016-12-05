library(OncoSimulR)



ng <- 50
u <- allFitnessEffects(noIntGenes = c(rep(0.1, ng/2), rep(-0.1, ng/2)))
t_e_10000 <- system.time(e_10000 <- oncoSimulPop(5,
                                                 u,
                                                 model = "Exp",
                                                 mu = 1e-7,
                                                 detectionSize = 1e9,
                                                 initSize = 1e5,
                                                 detectionDrivers = NA,
                                                 detectionProb = NA,
                                                 keepPhylog = TRUE,
                                                 onlyCancer = FALSE,
                                                 mutationPropGrowth = FALSE,
                                                 keepEvery = NA,
                                                 mc.cores = 1
                                ))
t_e_10000
summary(e_10000)[, c(1:3, 8, 9)]
print(object.size(e_10000), units = "MB")



## Work on this
ng <- 50
u <- allFitnessEffects(noIntGenes = c(rep(0.1, ng/2), rep(-0.1, ng/2)))
t_e_10000 <- system.time(e_10000 <- oncoSimulPop(5,
                                                 u,
                                                 model = "Exp",
                                                 mu = 1e-7,
                                                 detectionSize = 1e10,
                                                 initSize = 1e5,
                                                 detectionDrivers = NA,
                                                 detectionProb = NA,
                                                 keepPhylog = TRUE,
                                                 onlyCancer = FALSE,
                                                 mutationPropGrowth = FALSE,
                                                 keepEvery = NA,
                                                 mc.cores = 1
                                ))
t_e_10000
summary(e_10000)[, c(1:3, 8, 9)]
print(object.size(e_10000), units = "MB")











rm(list = ls())
gc()
ng <- 10
u <- allFitnessEffects(noIntGenes = c(rep(0.1, ng)))

## Avoid extinction so have all genes positive fitnes plus
## large init Size.


rm(list = ls())
gc()

ng <- 10
u <- allFitnessEffects(noIntGenes = c(rep(0.1, ng)))


t_e_7 <- system.time(e_7 <- oncoSimulPop(5,
                                         u,
                                         model = "Exp",
                                         mu = 1e-7,
                                         detectionSize = 1e7,
                                         initSize = 1e5,
                                         detectionDrivers = NA,
                                         detectionProb = NA,
                                         keepPhylog = TRUE,
                                         onlyCancer = FALSE,
                                         finalTime = 5000,
                                         mutationPropGrowth = FALSE,
                                         keepEvery = NA,
                                         mc.cores = 1
                                         ))
t_e_7
summary(e_7)[, c(1:3, 8, 9)]
print(object.size(e_7), units = "MB")


rm(list = ls())
gc()

ng <- 10
u <- allFitnessEffects(noIntGenes = c(rep(0.1, ng)))

t_e_8 <- system.time(e_8 <- oncoSimulPop(5,
                                         u,
                                         model = "Exp",
                                         mu = 1e-7,
                                         detectionSize = 1e8,
                                         initSize = 1e5,
                                         detectionDrivers = NA,
                                         detectionProb = NA,
                                         keepPhylog = TRUE,
                                         onlyCancer = FALSE,
                                         finalTime = 5000,
                                         mutationPropGrowth = FALSE,
                                         keepEvery = NA,
                                         mc.cores = 1
                                         ))
t_e_8
summary(e_8)[, c(1:3, 8, 9)]
print(object.size(e_8), units = "MB")


rm(list = ls())
gc()

ng <- 10
u <- allFitnessEffects(noIntGenes = c(rep(0.1, ng)))

t_e_9 <- system.time(e_9 <- oncoSimulPop(5,
                                         u,
                                         model = "Exp",
                                         mu = 1e-7,
                                         detectionSize = 1e9,
                                         initSize = 1e5,
                                         detectionDrivers = NA,
                                         detectionProb = NA,
                                         keepPhylog = TRUE,
                                         onlyCancer = FALSE,
                                         finalTime = 5000,
                                         mutationPropGrowth = FALSE,
                                         keepEvery = NA,
                                         mc.cores = 1
                                         ))
t_e_9
summary(e_9)[, c(1:3, 8, 9)]
print(object.size(e_9), units = "MB")


rm(list = ls())
gc()

ng <- 10
u <- allFitnessEffects(noIntGenes = c(rep(0.1, ng)))

t_e_10 <- system.time(e_10 <- oncoSimulPop(5,
                                         u,
                                         model = "Exp",
                                         mu = 1e-7,
                                         detectionSize = 1e10,
                                         initSize = 1e5,
                                         detectionDrivers = NA,
                                         detectionProb = NA,
                                         keepPhylog = TRUE,
                                         onlyCancer = FALSE,
                                         finalTime = 5000,
                                         mutationPropGrowth = FALSE,
                                         keepEvery = NA,
                                         mc.cores = 1
                                         ))
t_e_10
summary(e_10)[, c(1:3, 8, 9)]
print(object.size(e_10), units = "MB")


rm(list = ls())
gc()

ng <- 10
u <- allFitnessEffects(noIntGenes = c(rep(0.1, ng)))

t_e_105 <- system.time(e_105 <- oncoSimulPop(5,
                                         u,
                                         model = "Exp",
                                         mu = 1e-7,
                                         detectionSize = 5e10,
                                         initSize = 1e5,
                                         detectionDrivers = NA,
                                         detectionProb = NA,
                                         keepPhylog = TRUE,
                                         onlyCancer = FALSE,
                                         finalTime = 5000,
                                         mutationPropGrowth = FALSE,
                                         keepEvery = NA,
                                         mc.cores = 1
                                         ))
t_e_105
summary(e_105)[, c(1:3, 8, 9)]
print(object.size(e_105), units = "MB")


rm(list = ls())
gc()

ng <- 10
u <- allFitnessEffects(noIntGenes = c(rep(0.1, ng)))

t_e_11 <- system.time(e_11 <- oncoSimulPop(5,
                                         u,
                                         model = "Exp",
                                         mu = 1e-7,
                                         detectionSize = 1e11,
                                         initSize = 1e5,
                                         detectionDrivers = NA,
                                         detectionProb = NA,
                                         keepPhylog = TRUE,
                                         onlyCancer = FALSE,
                                         finalTime = 5000,
                                         mutationPropGrowth = FALSE,
                                         keepEvery = NA,
                                         mc.cores = 1
                                         ))
t_e_11
summary(e_11)[, c(1:3, 8, 9)]
print(object.size(e_11), units = "MB")

rm(list = ls())
gc()

ng <- 10
u <- allFitnessEffects(noIntGenes = c(rep(0.1, ng)))

t_e_112 <- system.time(e_112 <- oncoSimulPop(5,
                                         u,
                                         model = "Exp",
                                         mu = 1e-7,
                                         detectionSize = 2e11,
                                         initSize = 1e5,
                                         detectionDrivers = NA,
                                         detectionProb = NA,
                                         keepPhylog = TRUE,
                                         onlyCancer = FALSE,
                                         finalTime = 5000,
                                         mutationPropGrowth = FALSE,
                                         keepEvery = NA,
                                         mc.cores = 1
                                         ))
t_e_112
summary(e_112)[, c(1:3, 8, 9)]
print(object.size(e_112), units = "MB")






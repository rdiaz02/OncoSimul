rm(list = ls())
set.seed(NULL)

library(OncoSimulR)




### keepEvery = 1

rm(list = ls())

ng <- 50 
u <- allFitnessEffects(noIntGenes = c(rep(0.1, ng)))

t_mc_k_50_1e8 <- system.time(mc_k_50_1e8 <- oncoSimulPop(5,
                                                     u,
                                                     model = "McFL",
                                                     mu = 1e-7,
                                                     detectionSize = 1e8,
                                                     initSize = 1e5,
                                                     detectionDrivers = NA,
                                                     detectionProb = NA,
                                                     keepPhylog = TRUE,
                                                     onlyCancer = FALSE,
                                                     mutationPropGrowth = FALSE,
                                                     keepEvery = 1,
                                                     finalTime = 5000,
                                                     mc.cores = 1
                                                     ))
t_mc_k_50_1e8
try(summary(mc_k_50_1e8)[, c(1:3, 8, 9)])
print(object.size(mc_k_50_1e8), units = "MB")


t_mc_k_50_1e9 <- system.time(mc_k_50_1e9 <- oncoSimulPop(5,
                                                     u,
                                                     model = "McFL",
                                                     mu = 1e-7,
                                                     detectionSize = 1e9,
                                                     initSize = 1e5,
                                                     detectionDrivers = NA,
                                                     detectionProb = NA,
                                                     keepPhylog = TRUE,
                                                     onlyCancer = FALSE,
                                                     mutationPropGrowth = FALSE,
                                                     keepEvery = 1,
                                                     finalTime = 5000,
                                                     mc.cores = 1
                                                     ))
t_mc_k_50_1e9
try(summary(mc_k_50_1e9)[, c(1:3, 8, 9)])
print(object.size(mc_k_50_1e9), units = "MB")


t_mc_k_50_1e10 <- system.time(mc_k_50_1e10 <- oncoSimulPop(5,
                                                     u,
                                                     model = "McFL",
                                                     mu = 1e-7,
                                                     detectionSize = 1e10,
                                                     initSize = 1e5,
                                                     detectionDrivers = NA,
                                                     detectionProb = NA,
                                                     keepPhylog = TRUE,
                                                     onlyCancer = FALSE,
                                                     mutationPropGrowth = FALSE,
                                                     keepEvery = 1,
                                                     finalTime = 5000,
                                                     mc.cores = 1
                                                     ))
t_mc_k_50_1e10
try(summary(mc_k_50_1e10)[, c(1:3, 8, 9)])
print(object.size(mc_k_50_1e10), units = "MB")


t_mc_k_50_5e10 <- system.time(mc_k_50_5e10 <- oncoSimulPop(5,
                                                     u,
                                                     model = "McFL",
                                                     mu = 1e-7,
                                                     detectionSize = 5e10,
                                                     initSize = 1e5,
                                                     detectionDrivers = NA,
                                                     detectionProb = NA,
                                                     keepPhylog = TRUE,
                                                     onlyCancer = FALSE,
                                                     mutationPropGrowth = FALSE,
                                                     keepEvery = 1,
                                                     finalTime = 5000,
                                                     mc.cores = 1,
                                                     max.wall.time = 600
                                                     ))
t_mc_k_50_5e10
try(summary(mc_k_50_5e10)[, c(1:3, 8, 9)])
print(object.size(mc_k_50_5e10), units = "MB")


t_mc_k_50_1e11 <- system.time(mc_k_50_1e11 <- oncoSimulPop(5,
                                                     u,
                                                     model = "McFL",
                                                     mu = 1e-7,
                                                     detectionSize = 1e11,
                                                     initSize = 1e5,
                                                     detectionDrivers = NA,
                                                     detectionProb = NA,
                                                     keepPhylog = TRUE,
                                                     onlyCancer = FALSE,
                                                     mutationPropGrowth = FALSE,
                                                     keepEvery = 1,
                                                     finalTime = 5000,
                                                     mc.cores = 1,
                                                     max.wall.time = 600
                                                     ))
try(summary(mc_k_50_1e11)[, c(1:3, 8, 9)])
print(object.size(mc_k_50_1e11), units = "MB")


t_mc_k_50_5e11 <- system.time(mc_k_50_5e11 <- oncoSimulPop(5,
                                                     u,
                                                     model = "McFL",
                                                     mu = 1e-7,
                                                     detectionSize = 5e11,
                                                     initSize = 1e5,
                                                     detectionDrivers = NA,
                                                     detectionProb = NA,
                                                     keepPhylog = TRUE,
                                                     onlyCancer = FALSE,
                                                     mutationPropGrowth = FALSE,
                                                     keepEvery = 1,
                                                     finalTime = 5000,
                                                     mc.cores = 1,
                                                     max.wall.time = 600
                                                     ))
try(summary(mc_k_50_5e11)[, c(1:3, 8, 9)])
print(object.size(mc_k_50_5e11), units = "MB")




###########################################################



ng <- 50 
u <- allFitnessEffects(noIntGenes = c(rep(0.1, ng)))

t_exp_k_50_1e8 <- system.time(exp_k_50_1e8 <- oncoSimulPop(5,
                                                     u,
                                                     model = "Exp",
                                                     mu = 1e-7,
                                                     detectionSize = 1e8,
                                                     initSize = 1e5,
                                                     detectionDrivers = NA,
                                                     detectionProb = NA,
                                                     keepPhylog = TRUE,
                                                     onlyCancer = FALSE,
                                                     mutationPropGrowth = FALSE,
                                                     keepEvery = 1,
                                                     finalTime = 5000,
                                                     mc.cores = 1
                                                     ))
t_exp_k_50_1e8
try(summary(exp_k_50_1e8)[, c(1:3, 8, 9)])
print(object.size(exp_k_50_1e8), units = "MB")


t_exp_k_50_1e9 <- system.time(exp_k_50_1e9 <- oncoSimulPop(5,
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
                                                     keepEvery = 1,
                                                     finalTime = 5000,
                                                     mc.cores = 1
                                                     ))
t_exp_k_50_1e9
try(summary(exp_k_50_1e9)[, c(1:3, 8, 9)])
print(object.size(exp_k_50_1e9), units = "MB")


t_exp_k_50_1e10 <- system.time(exp_k_50_1e10 <- oncoSimulPop(5,
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
                                                     keepEvery = 1,
                                                     finalTime = 5000,
                                                     mc.cores = 1,
                                                     max.wall.time = 600
                                                     ))
t_exp_k_50_1e10
try(summary(exp_k_50_1e10)[, c(1:3, 8, 9)])
print(object.size(exp_k_50_1e10), units = "MB")



t_exp_k_50_5e10 <- system.time(exp_k_50_5e10 <- oncoSimulPop(5,
                                                     u,
                                                     model = "Exp",
                                                     mu = 1e-7,
                                                     detectionSize = 5e10,
                                                     initSize = 1e5,
                                                     detectionDrivers = NA,
                                                     detectionProb = NA,
                                                     keepPhylog = TRUE,
                                                     onlyCancer = FALSE,
                                                     mutationPropGrowth = FALSE,
                                                     keepEvery = 1,
                                                     finalTime = 5000,
                                                     mc.cores = 1,
                                                     max.wall.time = 600,
                                                     errorHitWallTime = FALSE,
                                                     errorHitMaxTries = FALSE
                                                     ))
t_exp_k_50_5e10
try(summary(exp_k_50_5e10)[, c(1:3, 8, 9)])
print(object.size(exp_k_50_5e10), units = "MB")



t_exp_k_50_1e11 <- system.time(exp_k_50_1e11 <- oncoSimulPop(5,
                                                     u,
                                                     model = "Exp",
                                                     mu = 1e-7,
                                                     detectionSize = 1e11,
                                                     initSize = 1e5,
                                                     detectionDrivers = NA,
                                                     detectionProb = NA,
                                                     keepPhylog = TRUE,
                                                     onlyCancer = FALSE,
                                                     mutationPropGrowth = FALSE,
                                                     keepEvery = 1,
                                                     finalTime = 5000,
                                                     mc.cores = 1,
                                                     max.wall.time = 600,
                                                     errorHitWallTime = FALSE,
                                                     errorHitMaxTries = FALSE
                                                     ))
t_exp_k_50_1e11
try(summary(exp_k_50_1e11)[, c(1:3, 8, 9)])
print(object.size(exp_k_50_1e11), units = "MB")


## Example code used in the vignette, but not executed there.

## This can use memory rather quickly. You might want to rm objects and
## gc().

library(OncoSimulR)
rm(list = ls()); gc()

ng <- 10000
u <- allFitnessEffects(noIntGenes = c(rep(0.1, ng/2), rep(-0.1, ng/2)))

t_e_10000 <- system.time(e_10000 <- oncoSimulPop(5,
                                                 u,
                                                 model = "Exp",
                                                 mu = 1e-7,
                                                 detectionSize = 1e6,
                                                 detectionDrivers = NA,
                                                 detectionProb = NA,
                                                 keepPhylog = TRUE,
                                                 onlyCancer = FALSE,
                                                 mutationPropGrowth = TRUE,
                                                 mc.cores = 1
                                ))
t_e_10000
summary(e_10000)[, c(1:3, 8, 9)]
print(object.size(e_10000), units = "MB")




t_e_10000b <- system.time(e_10000b <- oncoSimulPop(5,
                                                   u,
                                                   model = "Exp",
                                                   mu = 1e-7,
                                                   detectionSize = 1e6,
                                                   detectionDrivers = NA,
                                                   detectionProb = NA,
                                                   keepPhylog = TRUE,
                                                   onlyCancer = FALSE,
                                                   keepEvery = NA,
                                                   mutationPropGrowth = TRUE,
                                                   mc.cores = 1
                                ))
t_e_10000b
summary(e_10000b)[, c(1:3, 8, 9)]
print(object.size(e_10000b), units = "MB")



rm(list = ls()); gc()
ng <- 50000
u <- allFitnessEffects(noIntGenes = c(rep(0.1, ng/2), rep(-0.1, ng/2)))

t_e_50000 <- system.time(e_50000 <- oncoSimulPop(5,
                                                 u,
                                                 model = "Exp",
                                                 mu = 1e-7,
                                                 detectionSize = 1e6,
                                                 detectionDrivers = NA,
                                                 detectionProb = NA,
                                                 keepPhylog = TRUE,
                                                 onlyCancer = FALSE,
                                                 keepEvery = NA,
                                                 mutationPropGrowth = FALSE,
                                                 mc.cores = 1
                                                 ))

t_e_50000

summary(e_50000)[, c(1:3, 8, 9)]

print(object.size(e_50000), units = "MB")


rm(list = ls()); gc()
ng <- 50000
u <- allFitnessEffects(noIntGenes = c(rep(0.1, ng/2), rep(-0.1, ng/2)))

t_e_50000np <- system.time(e_50000np <- oncoSimulPop(5,
                                                 u,
                                                 model = "Exp",
                                                 mu = 1e-7,
                                                 detectionSize = 1e6,
                                                 detectionDrivers = NA,
                                                 detectionProb = NA,
                                                 keepPhylog = TRUE,
                                                 onlyCancer = FALSE,
                                                 keepEvery = 1,
                                                 mutationPropGrowth = FALSE,
                                                 mc.cores = 1
                                                 ))

t_e_50000np

summary(e_50000np)[, c(1:3, 8, 9)]

print(object.size(e_50000np), units = "MB")



rm(list = ls()); gc()
ng <- 50000
u <- allFitnessEffects(noIntGenes = c(rep(0.1, ng/2), rep(-0.1, ng/2)))

t_e_50000c <- system.time(e_50000c <- oncoSimulPop(5,
                                                 u,
                                                 model = "Exp",
                                                 mu = 1e-7,
                                                 detectionSize = 1e6,
                                                 detectionDrivers = NA,
                                                 detectionProb = NA,
                                                 keepPhylog = TRUE,
                                                 onlyCancer = FALSE,
                                                 keepEvery = NA,
                                                 mutationPropGrowth = TRUE,
                                                 mc.cores = 1
                                                 ))

t_e_50000c

summary(e_50000c)[, c(1:3, 8, 9)]

print(object.size(e_50000c), units = "MB")



#### McFL


rm(list = ls()); gc()
ng <- 50000
u <- allFitnessEffects(noIntGenes = c(rep(0.1, ng/2), rep(-0.1, ng/2)))

t_mc_50000_nmpg <- system.time(mc_50000_nmpg <- oncoSimulPop(5,
                                                   u,
                                                   model = "McFL",
                                                   mu = 1e-7,
                                                   detectionSize = 1e6,
                                                   detectionDrivers = NA,
                                                   detectionProb = NA,
                                                   keepPhylog = TRUE,
                                                   onlyCancer = FALSE,
                                                   keepEvery = NA,
                                                   mutationPropGrowth = FALSE,
                                                   mc.cores = 1
                                                   ))
t_mc_50000_nmpg

summary(mc_50000_nmpg)[, c(1:3, 8, 9)]

print(object.size(mc_50000_nmpg), units = "MB")



rm(list = ls()); gc()

ng <- 50000
u <- allFitnessEffects(noIntGenes = c(rep(0.1, ng/2), rep(-0.1, ng/2)))

t_mc_50000_nmpg_k <- system.time(mc_50000_nmpg_k <- oncoSimulPop(5,
                                                   u,
                                                   model = "McFL",
                                                   mu = 1e-7,
                                                   detectionSize = 1e6,
                                                   detectionDrivers = NA,
                                                   detectionProb = NA,
                                                   keepPhylog = TRUE,
                                                   onlyCancer = FALSE,
                                                   keepEvery = 1,
                                                   mutationPropGrowth = FALSE,
                                                   mc.cores = 1
                                                   ))
t_mc_50000_nmpg_k

summary(mc_50000_nmpg_k)[, c(1:3, 8, 9)]

print(object.size(mc_50000_nmpg_k), units = "MB")




rm(list = ls()); gc()
ng <- 50000
u <- allFitnessEffects(noIntGenes = c(rep(0.1, ng/2), rep(-0.1, ng/2)))

t_mc_50000_nmpg_3e6 <- system.time(mc_50000_nmpg_3e6 <- oncoSimulPop(5,
                                                   u,
                                                   model = "McFL",
                                                   mu = 1e-7,
                                                   detectionSize = 3e6,
                                                   detectionDrivers = NA,
                                                   detectionProb = NA,
                                                   keepPhylog = TRUE,
                                                   onlyCancer = FALSE,
                                                   keepEvery = NA,
                                                   mutationPropGrowth = FALSE,
                                                   mc.cores = 1
                                                   ))
t_mc_50000_nmpg_3e6

summary(mc_50000_nmpg_3e6)[, c(1:3, 8, 9)]

print(object.size(mc_50000_nmpg_3e6), units = "MB")



rm(list = ls()); gc()
ng <- 50000
u <- allFitnessEffects(noIntGenes = c(rep(0.1, ng/2), rep(-0.1, ng/2)))

t_mc_50000_nmpg_5mu <- system.time(mc_50000_nmpg_5mu <- oncoSimulPop(5,
                                                   u,
                                                   model = "McFL",
                                                   mu = 5e-7,
                                                   detectionSize = 1e6,
                                                   detectionDrivers = NA,
                                                   detectionProb = NA,
                                                   keepPhylog = TRUE,
                                                   onlyCancer = FALSE,
                                                   keepEvery = NA,
                                                   mutationPropGrowth = FALSE,
                                                   mc.cores = 1
                                                   ))
t_mc_50000_nmpg_5mu

summary(mc_50000_nmpg_5mu)[, c(1:3, 8, 9)]

print(object.size(mc_50000_nmpg_5mu), units = "MB")



rm(list = ls()); gc()
ng <- 50000
u <- allFitnessEffects(noIntGenes = c(rep(0.1, ng/2), rep(-0.1, ng/2)))

t_mc_50000_nmpg_5mu_k <- system.time(mc_50000_nmpg_5mu_k <- oncoSimulPop(5,
                                                   u,
                                                   model = "McFL",
                                                   mu = 5e-7,
                                                   detectionSize = 1e6,
                                                   detectionDrivers = NA,
                                                   detectionProb = NA,
                                                   keepPhylog = TRUE,
                                                   onlyCancer = FALSE,
                                                   keepEvery = 1,
                                                   mutationPropGrowth = FALSE,
                                                   mc.cores = 1
                                                   ))
t_mc_50000_nmpg_5mu_k

summary(mc_50000_nmpg_5mu_k)[, c(1:3, 8, 9)]

print(object.size(mc_50000_nmpg_5mu_k), units = "MB")





rm(list = ls()); gc()
ng <- 50000
u <- allFitnessEffects(noIntGenes = c(rep(0.1, ng/2), rep(-0.1, ng/2)))

t_mc_50000 <- system.time(mc_50000 <- oncoSimulPop(5,
                                                   u,
                                                   model = "McFL",
                                                   mu = 1e-7,
                                                   detectionSize = 1e6,
                                                   detectionDrivers = NA,
                                                   detectionProb = NA,
                                                   keepPhylog = TRUE,
                                                   onlyCancer = FALSE,
                                                   keepEvery = NA,
                                                   mutationPropGrowth = TRUE,
                                                   mc.cores = 1
                                                   ))
t_mc_50000

summary(mc_50000)[, c(1:3, 8, 9)]

print(object.size(mc_50000), units = "MB")











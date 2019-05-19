options(stringsAsFactors = FALSE)

r <- data.frame(Genotype = c("WT", "cucu", "cococ", "cucu, cococ"), 
                 Fitness = c("1 + 1.5*f_",
                             "5 + 3*(f_A + f_B + f_A_B)",
                             "5 + 3*(f_A + f_B + f_A_B)",
                             "7 + 5*(f_A + f_B + f_A_B)"))

afe <- allFitnessEffects(genotFitness = r, 
                         frequencyDependentFitness = TRUE,
                         frequencyType = "rel",
                         spPopSizes = c(5000, 2500, 2500, 500))

evalAllGenotypes(afe)

## pass spPopSizes to evalAllGenotyes, not in allFitnessEffects



r <- data.frame(Genotype = c("WT", "cucu", "cococ"), 
                 Fitness = c("1 + 1.5*f_",
                             "5 + 3*(f_A + f_B)",
                             "6 + 3*(f_A + f_B)"))
afe <- allFitnessEffects(genotFitness = r, 
                         frequencyDependentFitness = TRUE,
                         frequencyType = "rel",
                         spPopSizes = c(5000, 2500, 2500))
evalAllGenotypes(afe)


r2 <- data.frame(Genotype = c("WT", "cucu", "cococ"), 
                 Fitness = c("1 + 1.5*f_",
                             "5 + 3*(f_1 + f_2)",
                             "6 + 3*(f_1 + f_2)"))
afe2 <- allFitnessEffects(genotFitness = r2, 
                         frequencyDependentFitness = TRUE,
                         frequencyType = "rel",
                         spPopSizes = c(5000, 2500, 2500))
evalAllGenotypes(afe2)





r2 <- data.frame(Genotype = c("WT", "cucu", "cococ"), 
                 Fitness = c("1",
                             "1.5 + 1*(f_2)",
                             "2 + 1*(f_1)"))

afe3 <- allFitnessEffects(genotFitness = r2, 
                         frequencyDependentFitness = TRUE,
                         frequencyType = "rel")

mt <- allMutatorEffects(epistasis = c("cucu" = 2,
                                      "cococ" = 4))


set.seed(1)
s1 <- oncoSimulIndiv(afe3, 
                     model = "McFL", 
                     onlyCancer = FALSE, 
                     finalTime = 100,
                     mu = 1e-5,
                     initSize = 5000, 
                     keepPhylog = FALSE,
                     seed = NULL, 
                     errorHitMaxTries = FALSE, 
                     errorHitWallTime = FALSE)
plot(s1, show = "genotypes")


set.seed(1)
s2 <- oncoSimulIndiv(afe3,
                     muEF = mt,
                     model = "McFL", 
                     onlyCancer = FALSE, 
                     finalTime = 100,
                     mu = 1e-5,
                     initSize = 5000, 
                     keepPhylog = FALSE,
                     seed = NULL, 
                     errorHitMaxTries = FALSE, 
                     errorHitWallTime = FALSE)
plot(s2, show = "genotypes")


################


r4 <- data.frame(Genotype = c("WT", "A", "B", "A, B"), 
                 Fitness = c("1",
                             "1.4 + 1*(f_2)",
                             "1.4 + 1*(f_1)",
                             "1.6 + f_1 + f_2"))
afe4 <- allFitnessEffects(genotFitness = r4, 
                         frequencyDependentFitness = TRUE,
                         frequencyType = "rel")


set.seed(1)
s1 <- oncoSimulIndiv(afe4, 
                     model = "McFL", 
                     onlyCancer = FALSE, 
                     finalTime = 50,
                     mu = 1e-4,
                     initSize = 5000, 
                     keepPhylog = TRUE,
                     seed = NULL, 
                     errorHitMaxTries = FALSE, 
                     errorHitWallTime = FALSE)
plot(s1, show = "genotypes")


mt <- allMutatorEffects(epistasis = c("A" = 0.1,
                                      "B" = 10))
set.seed(1)
s2 <- oncoSimulIndiv(afe4,
                     muEF = mt,
                     model = "McFL", 
                     onlyCancer = FALSE, 
                     finalTime = 50,
                     mu = 1e-4,
                     initSize = 5000, 
                     keepPhylog = TRUE,
                     seed = NULL, 
                     errorHitMaxTries = FALSE, 
                     errorHitWallTime = FALSE)
plot(s2, show = "genotypes")


plotClonePhylog(s1, keepEvents = TRUE)
plotClonePhylog(s2, keepEvents = TRUE)


##########3


gffd <- data.frame(Genotype = c("WT", "A", "B", "C", "A, B"), 
                   Fitness = c("1",
                               "1.2 + 0.5 * f_1_2",
                               "1.4 - 0.5 * f_3",
                               "2.6 + 0.7*(log(f_1 + f_2)) + f_1_2",
                               "1.2 + sqrt(f_1 + f_3 + f_2)"))
afd <- allFitnessEffects(genotFitness = gffd, 
                         frequencyDependentFitness = TRUE,
                         frequencyType = "rel")
set.seed(2)
sfd <- oncoSimulIndiv(afd, 
                     model = "McFL", 
                     onlyCancer = FALSE, 
                     finalTime = 50,
                     mu = 1e-4,
                     initSize = 10000, 
                     keepPhylog = FALSE,
                     seed = NULL, 
                     errorHitMaxTries = FALSE, 
                     errorHitWallTime = FALSE)
plot(sfd, show = "genotypes")


## 
gffd2 <- data.frame(Genotype = c("WT", "A", "B", "C"), 
                   Fitness = c("1",
                               "0.9 + 0.3 * f_2 + 0.3 * f_",
                               "0.9 + 0.3 * f_1 + 0.3 * f_",
                               "0.2 + 2 * (f_1 + f_2)"
                               ## "0.4 + 5 * f_1 + 5 * f_1"
                               ))
afd2 <- allFitnessEffects(genotFitness = gffd2, 
                         frequencyDependentFitness = TRUE,
                         frequencyType = "rel")

sfd2 <- oncoSimulIndiv(afd2, 
                     model = "McFL", 
                     onlyCancer = FALSE, 
                     finalTime = 100,
                     mu = 1e-4,
                     initSize = 5000, 
                     keepPhylog = FALSE,
                     seed = NULL, 
                     errorHitMaxTries = FALSE, 
                     errorHitWallTime = FALSE)
plot(sfd2, show = "genotypes")


## Can produce funny oscillations
gffd3 <- data.frame(Genotype = c("WT", "A", "B", "A, B"), 
                   Fitness = c("1",
                               "1 + 0.2 * (n_2 > 10)",
                               ".9 + 1 * (n_1 > 10)",
                               "0"# "0.02 + "
                               ## "0.4 + 5 * f_1 + 5 * f_1"
                               ))
afd3 <- allFitnessEffects(genotFitness = gffd3, 
                         frequencyDependentFitness = TRUE,
                         frequencyType = "abs")
set.seed(1)
sfd3 <- oncoSimulIndiv(afd3,
                       model = "McFL", 
                       onlyCancer = FALSE, 
                       finalTime = 200,
                       mu = 1e-4,
                       initSize = 5000, 
                       keepPhylog = FALSE,
                       seed = NULL, 
                       errorHitMaxTries = FALSE, 
                       errorHitWallTime = FALSE)
plot(sfd3, show = "genotypes", type = "line")
plot(sfd3, show = "genotypes")
sfd3

## we cannot get collapse, because death rate never larger than birth rate.
## Need to incorporate the changes from Diego's TFG



set.seed(15)
sfd3 <- oncoSi
mulIndiv(afd3,
                       model = "McFL", 
                       onlyCancer = FALSE, 
                       finalTime = 500,
                       mu = 5e-3,
                       initSize = 1000, 
                       keepPhylog = FALSE,
                       seed = NULL, 
                       errorHitMaxTries = FALSE, 
                       errorHitWallTime = FALSE)
plot(sfd3, show = "genotypes")
sfd3



### 
data(examplesFitnessEffects)
set.seed(1)
tmp <-  oncoSimulIndiv(examplesFitnessEffects[["o3"]],
                       model = "McFL", 
                       mu = 5e-5,
                       detectionSize = 1e8, 
                       detectionDrivers = 3,
                       sampleEvery = 0.025,
                       max.num.tries = 10,
                       keepEvery = 5,
                       initSize = 2000,
                       finalTime = 3000,
                       seed = NULL,
                       onlyCancer = FALSE,
                       keepPhylog = TRUE)
tmp


##
r4 <- rfitness(4)
r4[, 5] <- c(2, rep(2, 15))
af11 <- allFitnessEffects(genotFitness = r4)
## all fitness 1

tmp2 <- oncoSimulIndiv(af11,
                       model = "McFL", 
                       mu = 5e-5,
                       detectionSize = 1e8, 
                       detectionDrivers = 3,
                       sampleEvery = 0.025,
                       max.num.tries = 10,
                       keepEvery = 5,
                       initSize = 10000,
                       finalTime = 3000,
                       K = 100, ## 5820
                       seed = NULL,
                       onlyCancer = FALSE,
                       keepPhylog = TRUE)
tmp2



###############

## Define fitness of the different genotypes
gffd <- data.frame(Genotype = c("WT", "A", "B", "C", "A, B"), 
                   Fitness = c("1 + 1.5 * f_1_2",
                               "1.3 + 1.5 * f_1_2",
                               "1.4",
                               "1.1 + 0.7*((f_1 + f_2) > 0.3) + f_1_2",
                               "1.2 + sqrt(f_1 + f_3 + f_2) - 0.3 * (f_1_2 > 0.5)"),
                   stringsAsFactors = FALSE)
evalAllGenotypes(allFitnessEffects(genotFitness = gffd, 
                         frequencyDependentFitness = TRUE,
                         frequencyType = "rel",
                         spPopSizes = c(100, 20, 20, 30, 0)))

evalAllGenotypes(allFitnessEffects(genotFitness = gffd, 
                         frequencyDependentFitness = TRUE,
                         frequencyType = "rel",
                         spPopSizes = c(100, 30, 40, 0, 10)))

evalAllGenotypes(allFitnessEffects(genotFitness = gffd, 
                         frequencyDependentFitness = TRUE,
                         frequencyType = "rel",
                         spPopSizes = c(100, 30, 40, 0, 100)))

afd <- allFitnessEffects(genotFitness = gffd, 
                         frequencyDependentFitness = TRUE,
                         frequencyType = "rel")

set.seed(1) ## for reproducibility
sfd <- oncoSimulIndiv(afd, 
                     model = "McFL", 
                     onlyCancer = FALSE, 
                     finalTime = 100,
                     mu = 1e-4,
                     initSize = 5000, 
                     keepPhylog = FALSE,
                     seed = NULL, 
                     errorHitMaxTries = FALSE, 
                     errorHitWallTime = FALSE)
plot(sfd, show = "genotypes")




############


rar <- data.frame(Genotype = c("WT", "A", "B", "C"), 
                 Fitness = c("1",
                             "1.1 + .3*f_2",
                             "1.2 + .4*f_1",
                             "1.0 + .5 * (f_1 + f_2)"))
afear <- allFitnessEffects(genotFitness = rar, 
                         frequencyDependentFitness = TRUE,
                         frequencyType = "rel",
                         spPopSizes = c(100, 200, 300, 400))
evalAllGenotypes(afear)


rar2 <- data.frame(Genotype = c("WT", "A", "B", "C"), 
                 Fitness = c("1",
                             "1.1 + .3*(n_2/N)",
                             "1.2 + .4*(n_1/N)",
                             "1.0 + .5 * ((n_1 + n_2)/N)"))
afear2 <- allFitnessEffects(genotFitness = rar2, 
                         frequencyDependentFitness = TRUE,
                         frequencyType = "abs",
                         spPopSizes = c(100, 200, 300, 400))
evalAllGenotypes(afear2)


set.seed(1)
tmp1 <- oncoSimulIndiv(afear, 
                       model = "McFL", 
                       onlyCancer = FALSE, 
                       finalTime = 30,
                       mu = 1e-4,
                       initSize = 5000, 
                       keepPhylog = FALSE,
                       seed = NULL, 
                       errorHitMaxTries = FALSE, 
                       errorHitWallTime = FALSE)

set.seed(1)
tmp2 <- oncoSimulIndiv(afear2, 
                       model = "McFL", 
                       onlyCancer = FALSE, 
                       finalTime = 30,
                       mu = 1e-4,
                       initSize = 5000, 
                       keepPhylog = FALSE,
                       seed = NULL, 
                       errorHitMaxTries = FALSE, 
                       errorHitWallTime = FALSE)
identical(print(tmp1), print(tmp2))


rar3 <- data.frame(Genotype = c("WT", "A", "B", "C"), 
                 Fitness = c("1",
                             "1.1 + .3*(n_2/N)",
                             "1.2 + .4*(n_1/N)",
                             "1.0 + .5 * ( n_1 > 20)"))
afear3 <- allFitnessEffects(genotFitness = rar3, 
                         frequencyDependentFitness = TRUE,
                         frequencyType = "abs",
                         spPopSizes = c(100, 200, 300, 400))
evalAllGenotypes(afear3)


set.seed(1)
tmp3 <- oncoSimulIndiv(afear3, 
                       model = "McFL", 
                       onlyCancer = FALSE, 
                       finalTime = 100,
                       mu = 1e-4,
                       initSize = 5000, 
                       keepPhylog = FALSE,
                       seed = NULL, 
                       errorHitMaxTries = FALSE, 
                       errorHitWallTime = FALSE)
plot(tmp3, show = "genotypes")

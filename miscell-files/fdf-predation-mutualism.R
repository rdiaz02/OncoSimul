## General consumer-resource models. See Otto and Day, for example, setion 3.4

options(stringsAsFactors = FALSE)

fe_consumer_resource <- function(r1, r2, K1, K2, a_12, a_21, awt = 1e-5,
                                 gt = c("WT", "S1", "S2")) {
    data.frame(Genotype = gt,
               Fitness = c(
                  paste0("max(0, 1 - ", awt, " * N)"),
                  paste0("1 + ", r1,
                         " * ( 1 - (n_1 + ", a_12, " * n_2)/", K1,
                         ")"),
                  paste0("1 + ", r2,
                         " * ( 1 - (n_2 + ", a_21, " * n_1)/", K2,
                         ")")
                  ))
}
## Show expressions for fitness
fe_consumer_resource("r1", "r2", "K1", "K2", "a_12", "a_21", "awt")
## Predator-prey
fe_pred_prey <-
    allFitnessEffects(
        genotFitness =
            fe_consumer_resource(1.51, 1.41, 10000, 2000,
                                 0.3, -0.15,
                                 gt = c("WT","prey", "Predator")),
        frequencyDependentFitness = TRUE,
        frequencyType = "abs")
pred_prey <- oncoSimulIndiv(fe_pred_prey,
                            model = "Exp",
                            onlyCancer = FALSE, 
                            finalTime = 100,
                            mu = 1e-4,
                            initSize = 5e5, 
                            keepPhylog = TRUE,
                            seed = NULL, 
                            errorHitMaxTries = FALSE, 
                            errorHitWallTime = FALSE)
pred_prey




plot(pred_prey, show = "genotypes")

plot(pred_prey, show = "genotypes",
     xlim = c(80, 100))

plot(pred_prey, show = "genotypes", type = "line",
     xlim = c(80, 100), ylim = c(1500, 12000))


fe_pred_prey <-
    allFitnessEffects(
        genotFitness =
            fe_consumer_resource(1.51, 1.41, 100, 70,
                                 0.3, -0.15,
                                 gt = c("WT","prey", "Predator")),
        frequencyDependentFitness = TRUE,
        frequencyType = "abs")
pred_prey <- oncoSimulIndiv(fe_pred_prey,
                            model = "Exp",
                            onlyCancer = FALSE, 
                            finalTime = 1000,
                            mu = 1e-3,
                            initSize = 10000, 
                            keepPhylog = TRUE,
                            seed = NULL, 
                            errorHitMaxTries = FALSE, 
                            errorHitWallTime = FALSE)
plot(pred_prey, show = "genotypes")

plot(pred_prey, show = "genotypes",
     xlim = c(900, 1000))







## Commensalism
fe_commens <-
    allFitnessEffects(
        genotFitness =
            fe_consumer_resource(1.2, 1.3, 5000, 20000,
                                 0, -0.2,
                                 gt = c("WT","A", "Commensal")),
        frequencyDependentFitness = TRUE,
        frequencyType = "abs")

commens <- oncoSimulIndiv(fe_commens,
                            model = "Exp",
                            onlyCancer = FALSE, 
                            finalTime = 100,
                            mu = 1e-4,
                            initSize = 40000, 
                            keepPhylog = TRUE,
                            seed = NULL, 
                            errorHitMaxTries = FALSE, 
                            errorHitWallTime = FALSE)

plot(commens, show = "genotypes")

plot(commens, show = "genotypes",
     xlim = c(80, 100))

plot(commens, show = "genotypes", type = "line",
     xlim = c(80, 100), ylim = c(2000, 22000))







######
## check predators go extinct? Not with the simple models as these
## because even when prey = 0, no extincting as positive growth rates if
## below K.

fe_consumer_resource <- function(r1, r2, K1, K2, a_12, a_21, awt = 1e-5,
                                 gt = c("WT", "S1", "S2")) {
    data.frame(Genotype = gt,
               Fitness = c(
                  paste0("max(0.1, 1 - ", awt, " * (n_2 + n_1))"),
                  paste0("1 + ", r1,
                         " * ( 1 - (n_1 + ", a_12, " * n_2)/", K1,
                         ")"),
                  paste0("1 + ", r2,
                         " * ( 1 - (n_2 + ", a_21, " * n_1)/", K2,
                         ")")
                  ))
}

## Show expressions for birth rates
fe_consumer_resource("r1", "r2", "K1", "K2", "a_12", "a_21", "awt")


fe_consumer_resource(1.5, 1.4, 100, 50, 0.3, -0.15, awt = 0.5)




fe_pred_prey <-
    allFitnessEffects(
        genotFitness =
            fe_consumer_resource(1.5, 1.4, 100, 50,
                                 0.3, -0.15, awt = 0.5,
                                 gt = c("WT","prey", "Predator")),
        frequencyDependentFitness = TRUE,
        frequencyType = "abs")



evalAllGenotypes(allFitnessEffects(
    genotFitness =
        fe_consumer_resource(1.5, 1.4, 100, 50,
                             0.3, -0.15, awt = 0.5,
                             gt = c("WT","prey", "Predator")),
    frequencyDependentFitness = TRUE,
    frequencyType = "abs",
    spPopSizes = c(0, 30, 50)))




####
set.seed(1)
fe_pred_prey <-
    allFitnessEffects(
        genotFitness =
            fe_consumer_resource(1.5, 1.4, 10000, 2000,
                                 0.3, -0.15,
                                 gt = c("WT","prey", "Predator")),
        frequencyDependentFitness = TRUE,
        frequencyType = "abs")
pred_prey <- oncoSimulIndiv(fe_pred_prey,
                            model = "Exp",
                            onlyCancer = FALSE, 
                            finalTime = 100,
                            mu = 1e-4,
                            initSize = 40000, 
                            keepPhylog = TRUE,
                            seed = NULL, 
                            errorHitMaxTries = FALSE, 
                            errorHitWallTime = FALSE)
plot(pred_prey, show = "genotypes",
     xlim = c(80, 100))
plot(pred_prey, show = "genotypes", type = "line",
     xlim = c(80, 100), ylim = c(1500, 12000))














## Comptetition model
C_fe_LV <- function(r1, r2, K1, K2, a_12, a_21, awt = 1e-4,
                                 gt = c("WT", "S1", "S2")) {
    data.frame(Genotype = gt,
               Fitness = c(
                  paste0("max(0.1, 1 - ", awt, " * (n_2 + n_1))"),
                  paste0("1 + ", r1,
                         " * ( 1 - (n_1 + ", a_12, " * n_2)/", K1,
                         ")"),
                  paste0("1 + ", r2,
                         " * ( 1 - (n_2 + ", a_21, " * n_1)/", K2,
                         ")")
                  ))
}

fe_competition <-
    allFitnessEffects(
        genotFitness =
            C_fe_LV(1.5, 1.4, 10000, 4000, 0.6, 0.2,
                  gt = c("WT","S1", "S2")),
        frequencyDependentFitness = TRUE,
        frequencyType = "abs")

set.seed(1)

competition <- oncoSimulIndiv(fe_competition,
                            model = "Exp",
                            onlyCancer = FALSE, 
                            finalTime = 100,
                            mu = 1e-4,
                            initSize = 40000, 
                            keepPhylog = TRUE,
                            seed = NULL, 
                            errorHitMaxTries = FALSE, 
                            errorHitWallTime = FALSE)
plot(competition, show = "genotypes")
plot(competition, show = "genotypes",
     xlim = c(80, 100))




##### Pred prey

set.seed(1)

fe_pred_prey <-
    allFitnessEffects(
        genotFitness =
            C_fe_LV(1.5, 1.4, 10000, 4000, 0.6, -0.5,
                  gt = c("WT","prey", "Predator")),
        frequencyDependentFitness = TRUE,
        frequencyType = "abs")

pred_prey <- oncoSimulIndiv(fe_pred_prey,
                            model = "Exp",
                            onlyCancer = FALSE, 
                            finalTime = 100,
                            mu = 1e-4,
                            initSize = 40000, 
                            keepPhylog = TRUE,
                            seed = NULL, 
                            errorHitMaxTries = FALSE, 
                            errorHitWallTime = FALSE)
plot(pred_prey, show = "genotypes")
plot(pred_prey, show = "genotypes",
     xlim = c(50, 100))



## smaller K
fe_pred_prey <-
    allFitnessEffects(
        genotFitness =
            G_fe_LV(1.5, 1.4, 100, 40, 0.6, -0.5, awt = 1,
                  gt = c("WT","prey", "Predator")),
        frequencyDependentFitness = TRUE,
        frequencyType = "abs")
set.seed(1)
pred_prey <- oncoSimulIndiv(fe_pred_prey,
                            model = "Exp",
                            onlyCancer = FALSE, 
                            finalTime = 200,
                            mu = 1e-3,
                            initSize = 1000, 
                            keepPhylog = TRUE,
                            seed = NULL, 
                            errorHitMaxTries = FALSE, 
                            errorHitWallTime = FALSE)
plot(pred_prey, show = "genotypes")


plot(pred_prey, show = "genotypes",
     xlim = c(50, 100))





evalAllGenotypes(allFitnessEffects(
    genotFitness =
        C_fe_LV(1.5, 1.4, 100, 40,
              0.6, -0.5, awt = 0.1,
              gt = c("WT","prey", "Predator")),
    frequencyDependentFitness = TRUE,
    frequencyType = "abs",
    spPopSizes = c(0, 0, 20)))


evalAllGenotypes(allFitnessEffects(
    genotFitness =
        C_fe_LV(1.5, 1.4, 100, 40,
              0.6, -0.5, awt = 0.1,
              gt = c("WT","prey", "Predator")),
    frequencyDependentFitness = TRUE,
    frequencyType = "abs",
    spPopSizes = c(0, 0, 40)))



## Pred-prey, p. 76  of otto and day

C_fe_pred_prey <- function(r, a, c, e, d, awt = 0.1,
                           gt = c("WT", "prey", "Predator")) {
    data.frame(Genotype = gt,
               Fitness = c(
                   paste0("max(0.1, 1 - ", awt,
                          " * (n_2 + n_1))"),
                   paste0("1 + ", r, " - ", a,
                          " * ", c, " * n_2"),
                   paste0("1 + ", e, " * ", a,
                          " * ", c, " * n_1 - ", d)
               ))
}

C_fe_pred_prey("r", "a", "c", "e", "d")



evalAllGenotypes( allFitnessEffects(
        genotFitness =
            C_fe_pred_prey(r = 5, a = 1, c = 0.005,
                           e = 0.02, d = 0.5, awt = 0.0001,
                           gt = c("WT","prey", "Predator")),
    frequencyDependentFitness = TRUE,
    frequencyType = "abs",
    spPopSizes = c(1600,0, 1e7)))
    

fe_pred_prey2 <-
    allFitnessEffects(
        genotFitness =
            C_fe_pred_prey(r = 1.5, a = 1, c = 0.005,
                           e = 0.02, d = 0.4, awt = 0.001,
                           gt = c("WT","prey", "Predator")),
        frequencyDependentFitness = TRUE,
        frequencyType = "abs")
set.seed(2)
pred_prey <- oncoSimulIndiv(fe_pred_prey2,
                            model = "Exp",
                            sampleEvery = 0.01,
                            mu = 1e-3,
                            onlyCancer = FALSE, 
                            finalTime = 200,
                            initSize = 50000, 
                            keepPhylog = TRUE,
                            seed = NULL, 
                            errorHitMaxTries = FALSE, 
                            errorHitWallTime = FALSE)
plot(pred_prey, show = "genotypes")
pred_prey



fe_pred_prey2 <-
    allFitnessEffects(
        genotFitness =
            C_fe_pred_prey(r = .7, a = 1, c = 0.005,
                           e = 0.02, d = 0.4, awt = 0.001,
                           gt = c("WT","prey", "Predator")),
        frequencyDependentFitness = TRUE,
        frequencyType = "abs")
pred_prey <- oncoSimulIndiv(fe_pred_prey2,
                            model = "Exp",
                            sampleEvery = 0.01,
                            mu = 1e-3,
                            onlyCancer = FALSE, 
                            finalTime = 200,
                            initSize = 50000, 
                            keepPhylog = TRUE,
                            seed = NULL, 
                            errorHitMaxTries = FALSE, 
                            errorHitWallTime = FALSE)
plot(pred_prey, show = "genotypes")
pred_prey



fe_pred_prey2 <-
    allFitnessEffects(
        genotFitness =
            C_fe_pred_prey(r = .5, a = 1, c = 0.01,
                           e = 0.02, d = 0.1, awt = 0.001,
                           gt = c("WT","prey", "Predator")),
        frequencyDependentFitness = TRUE,
        frequencyType = "abs")

set.seed(1)
pred_prey <- oncoSimulIndiv(fe_pred_prey2,
                            model = "Exp",
                            sampleEvery = 0.01,
                            mu = 1e-3,
                            onlyCancer = FALSE, 
                            finalTime = 80,
                            initSize = 1e4, 
                            keepPhylog = TRUE,
                            seed = NULL, 
                            errorHitMaxTries = FALSE, 
                            errorHitWallTime = FALSE)
plot(pred_prey, show = "genotypes")
pred_prey









fe_pred_prey2 <-
    allFitnessEffects(
        genotFitness =
            C_fe_pred_prey(r = .7, a = 1, c = 0.005,
                           e = 0.02, d = 0.4, awt = 0.001,
                           gt = c("WT","prey", "Predator")),
        frequencyDependentFitness = TRUE,
        frequencyType = "abs")

set.seed(2)


pred_prey <- oncoSimulIndiv(fe_pred_prey2,
                            model = "Exp",
                            sampleEvery = 0.01,
                            mu = 1e-3,
                            onlyCancer = FALSE, 
                            finalTime = 80,
                            initSize = 1e4, 
                            keepPhylog = TRUE,
                            seed = NULL, 
                            errorHitMaxTries = FALSE, 
                            errorHitWallTime = FALSE)
plot(pred_prey, show = "genotypes")
pred_prey





